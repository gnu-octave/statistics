## Copyright (C) 2026 Avanish Salunke <avanishsalunke16@gmail.com>
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

classdef LinearModel
  ## -*- texinfo -*-
  ## @deftp {statistics} LinearModel
  ##
  ## Linear regression model
  ##
  ## The @code{LinearModel} class represents a least-squares (or, optionally,
  ## robust) linear regression fit of a response variable to one or more
  ## predictor variables.  A @code{LinearModel} object is returned by the
  ## @code{fitlm} function and holds everything about the fit in one place:
  ## the fitted coefficients, the data and specification used to produce
  ## them, and the diagnostics needed to assess the quality of the fit.
  ##
  ## The properties of a @code{LinearModel} object fall into four groups:
  ##
  ## @multitable @columnfractions 0.22 0.76
  ## @headitem Group @tab Properties
  ##
  ## @item Coefficient estimates @tab @code{Coefficients} (a table of
  ## estimates, standard errors, t-statistics, and p-values for each term),
  ## @code{CoefficientCovariance}, @code{CoefficientNames}, and the
  ## coefficient counts @code{NumCoefficients} and
  ## @code{NumEstimatedCoefficients}.
  ##
  ## @item Summary statistics of the fit @tab @code{DFE},
  ## @code{Fitted}, @code{Residuals} (raw, Pearson, Studentized, and
  ## standardized), @code{Diagnostics} (leverage, Cook's distance, and other
  ## per-observation influence measures), @code{MSE}, @code{RMSE},
  ## @code{Rsquared} (ordinary and adjusted), @code{SSE}, @code{SSR},
  ## @code{SST}, @code{LogLikelihood}, @code{ModelCriterion} (AIC, BIC, etc.),
  ## and @code{ModelFitVsNullModel} (the F-test of the fitted model against an
  ## intercept-only model).
  ##
  ## @item Fitting method information @tab @code{Robust}, which records
  ## the weighting function and tuning constant used when the model is fit by
  ## robust regression, and is empty for an ordinary least squares fit, and
  ## @code{Steps}, which records the stepwise fitting information whenever
  ## the model was fit using stepwise regression.  This implementation
  ## currently always returns @code{Steps} as an empty structure.
  ##
  ## @item Input data properties @tab @code{Formula},
  ## @code{NumObservations}, @code{NumPredictors}, @code{NumVariables},
  ## @code{ObservationInfo} (which observations were used, excluded, missing,
  ## or weighted), @code{ObservationNames}, @code{PredictorNames},
  ## @code{ResponseName}, @code{VariableInfo}, @code{VariableNames}, and
  ## @code{Variables}.
  ## @end multitable
  ##
  ## A categorical predictor expands to indicator columns, one per level bar
  ## the reference level, which the intercept carries.  When the model has no
  ## intercept, the @emph{first} categorical predictor is given an indicator
  ## for every one of its levels instead, so that its coefficients are the
  ## group means; any further categorical predictor stays reference coded,
  ## which keeps the design full rank.  This differs from MATLAB, which omits
  ## the reference level whether or not an intercept is present and so cannot
  ## fit the reference group at all -- for a three-level grouping variable
  ## @code{g}, MATLAB fits @code{y ~ g - 1} with two coefficients, predicts
  ## exactly 0 for every observation in the omitted group, and reports a
  ## negative @math{R^2}.  This implementation returns three coefficients, one
  ## per group.
  ##
  ## A @code{LinearModel} object supports categorical predictors, which are
  ## automatically encoded internally as indicator (dummy) variables,
  ## observation weights for a weighted least squares fit, excluding specific
  ## observations from the fit, and robust regression using iteratively
  ## reweighted least squares.  Once fitted, the following methods are
  ## available on a @code{LinearModel} object:
  ##
  ## @multitable @columnfractions 0.2 0.78
  ## @headitem Method @tab Description
  ##
  ## @item @code{predict} @tab Predict responses at new predictor values
  ## given in a matrix or table, or reproduce the training fitted values when
  ## called with no new data.  Can also return pointwise or simultaneous
  ## confidence or prediction intervals alongside the point predictions.
  ##
  ## @item @code{feval} @tab Predict responses given predictors as
  ## separate scalar or vector arguments (one per predictor variable) instead
  ## of a single matrix, so a @code{LinearModel} object can be evaluated the
  ## same way as a plain function handle.  Returns point predictions only.
  ##
  ## @item @code{random} @tab Simulate new response values at new
  ## predictor locations by adding independent Gaussian noise, drawn from the
  ## estimated error variance @code{MSE}, to the fitted response.
  ##
  ## @item @code{coefCI} @tab Return Wald confidence intervals for every
  ## fitted coefficient at a chosen significance level (default @math{0.05}).
  ##
  ## @item @code{coefTest} @tab Test a linear hypothesis on the fitted
  ## coefficients.  With no arguments, tests the overall model F-test that
  ## all non-intercept coefficients are zero; a custom hypothesis can be
  ## given as a contrast matrix and, if needed, right-hand-side values.
  ## Returns the p-value, and optionally the F-statistic and its numerator
  ## degrees of freedom.
  ##
  ## @item @code{dwtest} @tab Durbin-Watson test for first-order
  ## autocorrelation among the model residuals, with a choice of exact or
  ## approximate p-value computation and a one- or two-sided alternative.
  ##
  ## @item @code{addTerms} @tab Return a new, refitted @code{LinearModel}
  ## with terms added to the current model specification, given as a
  ## Wilkinson formula fragment or a terms matrix.  Weights, excluded rows,
  ## and categorical encodings carry over automatically; the original model
  ## object is left unmodified.
  ##
  ## @item @code{removeTerms} @tab Return a new, refitted
  ## @code{LinearModel} with terms removed from the current model
  ## specification, given as a Wilkinson formula fragment or a terms matrix.
  ## Weights, excluded rows, and categorical encodings carry over
  ## automatically; the original model object is left unmodified.
  ##
  ## @item @code{plotResiduals} @tab Plot the model residuals.  Default
  ## is a probability density histogram; other supported plot types are
  ## @qcode{'fitted'}, @qcode{'caseorder'}, @qcode{'lagged'},
  ## @qcode{'probability'}, and @qcode{'observed'}.
  ##
  ## @item @code{plotDiagnostics} @tab Plot per-observation influence
  ## diagnostics.  Default is leverage by observation row number; other
  ## supported plot types are @qcode{'cookd'}, @qcode{'covratio'},
  ## @qcode{'dfbetas'}, @qcode{'dffits'}, @qcode{'s2_i'}, and
  ## @qcode{'contour'} (standardized residuals against leverage with Cook's
  ## distance contours).
  ##
  ## @item @code{plotEffects} @tab Plot the estimated main effect and
  ## 95% confidence interval of each predictor, evaluated between its
  ## observed minimum and maximum with all other predictors held at their
  ## observed means.
  ##
  ## @item @code{plotAdjustedResponse} @tab Plot the fitted response
  ## against a single predictor, with the other predictors averaged out by
  ## averaging the fitted values over the observations used in the fit.
  ##
  ## @item @code{plotAdded} @tab Plot the incremental effect of one or
  ## more terms on the response, after removing the effects of all other
  ## terms, along with the fitted line and its 95% confidence bounds.
  ##
  ## @item @code{plot} @tab Plot a default view of the model.  Creates an
  ## added variable plot for the whole model when more than one predictor
  ## is included, a scatter plot of the data with a fitted curve and 95%
  ## confidence bounds when exactly one predictor is included, or a
  ## histogram of the residuals when no predictors are included.
  ##
  ## @item @code{plotInteraction} @tab Plot the main and conditional effects
  ## of two predictors, or the adjusted response as a function of one
  ## predictor for several fixed values of the other, to visualize whether
  ## the two predictors interact.
  ##
  ## @item @code{compact} @tab Return a @code{CompactLinearModel} that
  ## discards the training data and per-observation diagnostics while
  ## retaining the coefficient estimates and fit statistics needed for
  ## prediction and inference.
  ##
  ## @item @code{anova} @tab Analysis of variance for the fitted model,
  ## reporting either the per-term breakdown of sums of squares or a
  ## summary table of the model against the total and residual variation.
  ##
  ## @item @code{step} @tab Improve the fitted model by one or more
  ## steps of stepwise term selection, returning a new, refitted
  ## @code{LinearModel} without modifying the original.
  ## @end multitable
  ##
  ## Create a @code{LinearModel} object by using the @code{fitlm} function or
  ## the class constructor directly.
  ##
  ## @seealso{fitlm}
  ## @end deftp

  properties(GetAccess = public, SetAccess = protected)

    ## Coefficient estimate properties

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} CoefficientCovariance
    ##
    ## Covariance matrix of coefficient estimates
    ##
    ## A @math{p}-by-@math{p} numeric matrix of covariance values for the
    ## coefficient estimates, where @math{p} is the number of coefficients in
    ## the fitted model as given by @code{NumCoefficients}.  This property is
    ## read-only.
    ##
    ## @end deftp
    CoefficientCovariance = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} CoefficientNames
    ##
    ## Coefficient names
    ##
    ## A cell array of character vectors, each containing the name of the
    ## corresponding model term (e.g., @qcode{'(Intercept)'}, @qcode{'x1'},
    ## @qcode{'x1:x2'}).  This property is read-only.
    ##
    ## @end deftp
    CoefficientNames = {};

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} Coefficients
    ##
    ## Coefficient values
    ##
    ## A table with one row for each coefficient and four columns:
    ## @itemize
    ## @item @code{Estimate} - estimated coefficient value
    ## @item @code{SE} - standard error of the estimate
    ## @item @code{tStat} - t-statistic for a two-sided test
    ## @item @code{pValue} - p-value for the t-statistic
    ## @end itemize
    ## Coefficients that are dropped due to rank deficiency have
    ## @code{Estimate = 0}, @code{SE = 0}, @code{tStat = NaN},
    ## @code{pValue = NaN}.  This property is read-only.
    ##
    ## @end deftp
    Coefficients = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} NumCoefficients
    ##
    ## Number of model coefficients
    ##
    ## A positive integer giving the total number of coefficients in the fitted
    ## model, including any coefficients set to zero because the model terms are
    ## rank deficient.  This property is read-only.
    ##
    ## @end deftp
    NumCoefficients = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} NumEstimatedCoefficients
    ##
    ## Number of estimated coefficients
    ##
    ## A positive integer giving the number of coefficients actually estimated,
    ## i.e., not set to zero due to rank deficiency.
    ## @code{NumEstimatedCoefficients} equals the degrees of freedom for
    ## regression.  This property is read-only.
    ##
    ## @end deftp
    NumEstimatedCoefficients = [];

    ## Summary statistic properties

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} DFE
    ##
    ## Degrees of freedom for error
    ##
    ## A positive integer equal to the number of observations minus the number
    ## of estimated coefficients: @code{DFE = NumObservations -
    ## NumEstimatedCoefficients}.  This property is read-only.
    ##
    ## @end deftp
    DFE = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} Diagnostics
    ##
    ## Observation diagnostics
    ##
    ## A table with one row per observation and seven columns:
    ## @itemize
    ## @item @code{Leverage} - diagonal of the hat matrix @math{H}
    ## @item @code{CooksDistance} - Cook's distance, a measure of scaled
    ##   change in fitted values
    ## @item @code{Dffits} - delete-1 scaled differences in fitted values
    ## @item @code{S2_i} - delete-1 residual variance estimate
    ## @item @code{CovRatio} - ratio of the determinant of the coefficient
    ##   covariance matrix with and without each observation
    ## @item @code{Dfbetas} - @math{n}-by-@math{p} matrix of scaled changes
    ##   in coefficient estimates when each observation is deleted in turn
    ## @item @code{HatMatrix} - @math{n}-by-@math{n} projection matrix such
    ##   that @code{Fitted = HatMatrix * y}
    ## @end itemize
    ## Rows not used in fitting have @code{NaN} in @code{CooksDistance},
    ## @code{Dffits}, @code{S2_i}, and @code{CovRatio}, and zeros in
    ## @code{Leverage}, @code{Dfbetas}, and @code{HatMatrix}.  This property
    ## is read-only.
    ##
    ## @end deftp
    Diagnostics = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} Fitted
    ##
    ## Fitted response values
    ##
    ## An @math{n}-by-1 numeric vector of predicted response values based on
    ## the training data, where @math{n} is the total number of observations
    ## including excluded and missing rows, which contain @code{NaN}.  Use
    ## @code{predict} to obtain predictions for new data or to compute
    ## confidence bounds.  This property is read-only.
    ##
    ## @end deftp
    Fitted = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} LogLikelihood
    ##
    ## Log-likelihood of the fitted model
    ##
    ## A scalar numeric value equal to the log-likelihood of the response
    ## values, assuming each response is normally distributed with mean equal
    ## to the fitted value and variance equal to @math{SSE/n} (the MLE
    ## variance estimate).  This property is read-only.
    ##
    ## @end deftp
    LogLikelihood = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} ModelCriterion
    ##
    ## Model comparison criteria
    ##
    ## A structure with four fields:
    ## @itemize
    ## @item @code{AIC} - Akaike information criterion:
    ##   @math{-2 * logL + 2 * m}
    ## @item @code{AICc} - AIC corrected for sample size:
    ##   @math{AIC + (2*m*(m+1))/(n-m-1)}
    ## @item @code{BIC} - Bayesian information criterion:
    ##   @math{-2 * logL + m * log(n)}
    ## @item @code{CAIC} - Consistent AIC:
    ##   @math{-2 * logL + m * (log(n) + 1)}
    ## @end itemize
    ## Here @math{logL} is @code{LogLikelihood}, @math{m} is
    ## @code{NumEstimatedCoefficients}, and @math{n} is
    ## @code{NumObservations}.  This property is read-only.
    ##
    ## @end deftp
    ModelCriterion = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} ModelFitVsNullModel
    ##
    ## F-test of the fitted model versus the null model
    ##
    ## A structure with three fields:
    ## @itemize
    ## @item @code{Fstat} - F-statistic of the fitted model versus a null
    ##   model containing only a constant term
    ## @item @code{Pvalue} - p-value for the F-statistic
    ## @item @code{NullModel} - character vector describing the null model
    ## @end itemize
    ## This property is read-only.
    ##
    ## @end deftp
    ModelFitVsNullModel = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} MSE
    ##
    ## Mean squared error
    ##
    ## A scalar numeric value equal to @math{SSE / DFE}, where @code{SSE} is
    ## the sum of squared errors and @code{DFE} is the degrees of freedom for
    ## error.  This property is read-only.
    ##
    ## @end deftp
    MSE = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} Residuals
    ##
    ## Residuals for the fitted model
    ##
    ## A table with one row per observation and four columns:
    ## @itemize
    ## @item @code{Raw} - observed minus fitted values
    ## @item @code{Pearson} - raw residuals divided by @code{RMSE}
    ## @item @code{Standardized} - internally studentized residuals; raw
    ##   residuals divided by their estimated standard deviation using the
    ##   full-model @code{MSE}
    ## @item @code{Studentized} - externally studentized residuals; each raw
    ##   residual divided by an estimate of the standard deviation based on
    ##   all observations except that one, using the delete-1 @code{S2_i}
    ## @end itemize
    ## Rows not used in the fit contain @code{NaN}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Residuals = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} RMSE
    ##
    ## Root mean squared error
    ##
    ## A scalar numeric value equal to @math{sqrt(MSE)}.  This property is
    ## read-only.
    ##
    ## @end deftp
    RMSE = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} Rsquared
    ##
    ## R-squared goodness-of-fit statistics
    ##
    ## A structure with two fields:
    ## @itemize
    ## @item @code{Ordinary} - coefficient of determination:
    ##   @math{R^2 = SSR / SST}
    ## @item @code{Adjusted} - adjusted @math{R^2} that accounts for the
    ##   number of coefficients in the model
    ## @end itemize
    ## This property is read-only.
    ##
    ## @end deftp
    Rsquared = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} SSE
    ##
    ## Sum of squared errors
    ##
    ## A scalar numeric value equal to the sum of squared residuals.  For a
    ## model with an intercept, @math{SST = SSE + SSR}.  For weighted fits,
    ## this is the weighted sum of squares.  This property is read-only.
    ##
    ## @end deftp
    SSE = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} SSR
    ##
    ## Regression sum of squares
    ##
    ## A scalar numeric value equal to the sum of squared deviations of the
    ## fitted values from the mean of the response.  For a model with an
    ## intercept, @math{SST = SSE + SSR}.  For weighted fits, this is the
    ## weighted sum of squares.  This property is read-only.
    ##
    ## @end deftp
    SSR = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} SST
    ##
    ## Total sum of squares
    ##
    ## A scalar numeric value equal to the sum of squared deviations of the
    ## response from its mean. For a model with an intercept,
    ## @math{SST = SSE + SSR}. For a robust fit, @math{SST = SSE + SSR}
    ## rather than the deviation from the mean. For weighted fits, this is
    ## the weighted sum of squares. This property is read-only.
    ##
    ## @end deftp
    SST = [];


    ## Fitting method properties

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} Robust
    ##
    ## Robust fit information
    ##
    ## A structure with three fields:
    ## @itemize
    ## @item @code{WgtFun} - robust weighting function name, e.g.
    ##   @qcode{'bisquare'}
    ## @item @code{Tune} - tuning constant; empty if @code{WgtFun} is
    ##   @qcode{'ols'} or a function handle with the default tuning constant
    ## @item @code{Weights} - vector of final iteration weights; empty for
    ##   a @code{CompactLinearModel} object
    ## @end itemize
    ## This structure is empty unless the model was fit using robust
    ## regression.  This property is read-only.
    ##
    ## @end deftp
    Robust = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} Steps
    ##
    ## Stepwise fitting information
    ##
    ## In MATLAB, this is a structure with seven fields (@code{Start},
    ## @code{Lower}, @code{Upper}, @code{Criterion}, @code{PEnter},
    ## @code{PRemove}, and @code{History}, the last being a table with one
    ## row per step and columns @code{Action}, @code{TermName}, @code{Terms},
    ## @code{DF}, @code{delDF}, @code{FStat}, @code{PValue}), populated
    ## whenever the model was fit using stepwise regression, and empty
    ## otherwise.
    ##
    ## This implementation always returns an empty struct, regardless of
    ## whether the model was fit with @code{stepwiselm} or @code{step}; the
    ## per-step history is not yet tracked.  This is a known deviation from
    ## MATLAB and not a deliberate design choice.  This property is
    ## read-only.
    ##
    ## @end deftp
    Steps = struct ();


    ## Input data properties

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} Formula
    ##
    ## Model formula information
    ##
    ## A @code{LinearFormula} object representing the model formula, with
    ## properties including @code{ResponseName}, @code{LinearPredictor},
    ## @code{PredictorNames}, @code{TermNames}, @code{HasIntercept},
    ## @code{Terms} (the terms matrix), and @code{InModel}.  Converting it with
    ## @code{char} renders the whole formula.  This property is read-only.
    ##
    ## @end deftp
    Formula = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} NumObservations
    ##
    ## Number of observations used in the fit
    ##
    ## A positive integer giving the number of observations actually used in
    ## fitting.  Rows with missing values and rows excluded via the
    ## @code{'Exclude'} name-value argument are not counted.  This property
    ## is read-only.
    ##
    ## @end deftp
    NumObservations = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} NumPredictors
    ##
    ## Number of predictor variables
    ##
    ## A positive integer giving the number of predictor variables used to
    ## fit the model.  This property is read-only.
    ##
    ## @end deftp
    NumPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} NumVariables
    ##
    ## Number of variables in the input data
    ##
    ## A positive integer giving the total number of variables in the input
    ## data, counting predictors, the response, and any unused columns.
    ## This property is read-only.
    ##
    ## @end deftp
    NumVariables = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} ObservationInfo
    ##
    ## Per-observation metadata
    ##
    ## An @math{n}-by-4 table where @math{n} is the total number of rows in
    ## the input data.  The four columns are:
    ## @itemize
    ## @item @code{Weights} - observation weight, default is 1
    ## @item @code{Excluded} - logical; true if excluded via the
    ##   @code{'Exclude'} argument
    ## @item @code{Missing} - logical; true if the row contains any
    ##   @code{NaN} value
    ## @item @code{Subset} - logical; true if the observation was used in
    ##   the fit, i.e. not excluded and not missing
    ## @end itemize
    ## This property is read-only.
    ##
    ## @end deftp
    ObservationInfo = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} ObservationNames
    ##
    ## Observation names
    ##
    ## A cell array of character vectors containing the names of the
    ## observations.  If the fit was based on a table that has row names,
    ## this property holds those names.  Otherwise it is an empty cell array.
    ## This property is read-only.
    ##
    ## @end deftp
    ObservationNames = {};

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} PredictorNames
    ##
    ## Names of predictor variables
    ##
    ## A cell array of character vectors containing the names of the
    ## predictor variables used to fit the model.  This property is
    ## read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector containing the name of the response variable.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName = '';

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} VariableInfo
    ##
    ## Information about input variables
    ##
    ## A table with one row per variable including any unused variables, and
    ## four columns:
    ## @itemize
    ## @item @code{Class} - variable class as a character vector, e.g.
    ##   @qcode{'double'} or @qcode{'categorical'}
    ## @item @code{Range} - for continuous variables, a two-element vector
    ##   @code{[min, max]}; for categorical variables, a vector of the
    ##   distinct values
    ## @item @code{InModel} - logical; true if the variable is in the
    ##   fitted model
    ## @item @code{IsCategorical} - logical; true if the variable is
    ##   categorical
    ## @end itemize
    ## This property is read-only.
    ##
    ## @end deftp
    VariableInfo = [];

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} VariableNames
    ##
    ## Names of all variables in the input data
    ##
    ## A cell array of character vectors containing the names of all
    ## variables, including predictors, the response, and unused variables.
    ## For table input these are the table column names.  For matrix input
    ## these are the values given by @code{'VarNames'}, defaulting to
    ## @qcode{@{'x1','x2',...,'xp','y'@}}.  This property is read-only.
    ##
    ## @end deftp
    VariableNames = {};

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} Variables
    ##
    ## Input data as a table
    ##
    ## A table containing predictor and response values for all observations,
    ## including unused variables.  For table input this is the full input
    ## table.  For matrix input this is a table constructed from the
    ## predictor matrix and response vector.  This property is read-only.
    ##
    ## @end deftp
    Variables = [];

  endproperties

  properties(Access = private, Hidden)

    ## Full design matrix, n by p_design, used for predictions
    DesignMatrix = [];

    ## Column indices of active coefficients in the design matrix
    ActiveCols = [];

    ## Whether the model includes an intercept term
    HasIntercept = true;

    ## Response vector, full n by 1 with NaN for non-subset rows
    ResponseVector = [];

    ## Full n by 1 observation weights
    WeightVector = [];

    ## n by 1 logical mask: true for rows used in the fit
    SubsetMask = [];

    ## Terms matrix from modelspec or parse_modelspec
    TermsMatrix = [];

    ## Categorical level info for re-encoding in predict
    CatLevelInfo = [];

    ## Predictor names after categorical dummy expansion
    EncPredictorNames = {};

    ## Cached per-predictor design contrasts used by plotEffects
    EffectContrasts = [];

    ## Cached per-predictor-pair design contrasts used by plotInteraction
    InteractionContrasts = [];

    ## Cached conceptual-term name and design-column grouping used by anova
    TermGroups = {};

    ## Encoded predictor matrix (Path B only), cached for refit
    EncodedPredMatrix = [];

    ## Parsed NV options stored for refit operations
    OrigOpts = [];

  endproperties

  methods(Hidden)

    ## Custom display
    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ("%s =\n", in_name);
      endif
      disp (this);
    endfunction

    ## Custom display
    function disp (this)
      if (isempty (this.Robust))
        fprintf ("\n  Linear regression model:\n");
      else
        fprintf ("\n  Linear regression model (robust fit):\n");
      endif
      if (! isempty (this.Formula) && isa (this.Formula, 'LinearFormula'))
        fprintf ("      %s\n", char (this.Formula));
      endif

      if (! isempty (this.Coefficients))
        fprintf ("\n  Estimated Coefficients:\n\n");
        disp (this.Coefficients);
      endif

      fprintf ("\n");
      if (! isempty (this.NumObservations) && ! isempty (this.DFE))
        fprintf ("Number of observations: %d, Error degrees of freedom: %d\n", ...
                 this.NumObservations, this.DFE);
      endif
      if (! isempty (this.RMSE))
        fprintf ("Root Mean Squared Error: %g\n", this.RMSE);
      endif
      if (! isempty (this.Rsquared) && isstruct (this.Rsquared))
        fprintf ("R-squared: %g,  Adjusted R-Squared: %g\n", ...
                 this.Rsquared.Ordinary, this.Rsquared.Adjusted);
      endif
      if (! isempty (this.ModelFitVsNullModel) ...
          && isstruct (this.ModelFitVsNullModel) ...
          && isfield (this.ModelFitVsNullModel, 'Fstat'))
        fprintf ("F-statistic vs. constant model: %g, p-value = %g\n", ...
                 this.ModelFitVsNullModel.Fstat, ...
                 this.ModelFitVsNullModel.Pvalue);
      endif
    endfunction

    ## Class specific subscripted reference
    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case '()'
          error (strcat ("LinearModel: () indexing is not supported.", "  Use dot notation to access properties."));
        case '{}'
          error (strcat ("LinearModel: {} indexing is not supported.", "  Use dot notation to access properties."));
        case '.'
          if (! ischar (s.subs))
            error ("LinearModel.subsref: property name must be a character vector.");
          endif

          ## Allow normal execution if the user is calling a class method
          if (ismethod (this, s.subs))
            [varargout{1:nargout}] = builtin ('subsref', this, [s, chain_s]);
            return;
          endif
          try
            out = this.(s.subs);
          catch
            error ("LinearModel.subsref: unknown property '%s'.", s.subs);
          end_try_catch
      endswitch
      if (! isempty (chain_s))
        out = subsref (out, chain_s);
      endif
      varargout{1} = out;
    endfunction

  endmethods

  methods(Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {@var{mdl} =} LinearModel (@var{X}, @var{y})
    ## @deftypefnx {LinearModel} {@var{mdl} =} LinearModel (@var{tbl}, @var{resp_input})
    ## @deftypefnx {LinearModel} {@var{mdl} =} LinearModel (@dots{}, @var{modelspec})
    ## @deftypefnx {LinearModel} {@var{mdl} =} LinearModel (@dots{}, @var{Name}, @var{Value}, @dots{})
    ##
    ## Create a @qcode{LinearModel} class object representing a linear
    ## regression model.
    ##
    ## @code{@var{mdl} = LinearModel (@var{X}, @var{y})} returns a
    ## @code{LinearModel} object fit to the response @var{y} and the
    ## predictor data @var{X}.  Unless removed via the @qcode{'Intercept'}
    ## option, the fitted model contains a constant (intercept) term and one
    ## linear term for every column of @var{X}.
    ##
    ## @itemize
    ## @item
    ## @var{X} is an @math{N*P} numeric or logical matrix of predictor data,
    ## where rows correspond to observations and columns correspond to
    ## variables.  By default, the predictors are named @qcode{'x1'},
    ## @qcode{'x2'}, @dots{}, @qcode{'xP'}.
    ## @item
    ## @var{y} is an @math{N*1} numeric or logical vector of response values,
    ## and must have the same number of observations (rows) as @var{X}.  By
    ## default, the response is named @qcode{'y'}.
    ## @end itemize
    ##
    ## @code{@var{mdl} = LinearModel (@var{tbl}, @var{resp_input})} fits a
    ## model using the variables in the table (or dataset) @var{tbl} as
    ## predictors.  @var{resp_input} selects the response and can be a
    ## character vector naming a variable in @var{tbl}, or a numeric vector
    ## the same height as @var{tbl} to use as an external response.  If
    ## @var{resp_input} is left empty, the last variable in @var{tbl} is used
    ## as the response.  Variables that are @code{categorical} arrays, cell
    ## arrays of character vectors, or logical arrays are automatically
    ## treated as categorical predictors.
    ##
    ## @code{@var{mdl} = LinearModel (@dots{}, @var{modelspec})} additionally
    ## specifies the terms of the model to fit.  @var{modelspec} can be any of
    ## the following.
    ##
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Value} @tab @var{Description}
    ##
    ## @item @qcode{'constant'} @tab Model contains only an intercept
    ## term.
    ##
    ## @item @qcode{'linear'} @tab Model contains an intercept and one
    ## term for each predictor variable.  This is the default when
    ## @var{modelspec} is not specified.
    ##
    ## @item @qcode{'interactions'} @tab Model contains an intercept, all
    ## linear terms, and all pairwise products of distinct predictor
    ## variables (no squared terms).
    ##
    ## @item @qcode{'purequadratic'} @tab Model contains an intercept,
    ## all linear terms, and all squared terms.
    ##
    ## @item @qcode{'quadratic'} @tab Model contains an intercept, all
    ## linear terms, all pairwise products of distinct predictor variables,
    ## and all squared terms.
    ##
    ## @item @qcode{'full'} @tab Model contains an intercept and all
    ## terms up to and including the full @math{P}-way interaction of the
    ## predictor variables.
    ##
    ## @item terms matrix @tab A @math{T*P} or @math{T*(P+1)} numeric
    ## matrix, where @math{T} is the number of terms and @math{P} is the
    ## number of predictor variables.  Each row represents one term, and the
    ## value in column @math{j} is the exponent to which predictor @math{j}
    ## is raised in that term; a row of all zeros represents the intercept.
    ## If a @math{T*(P+1)} matrix is supplied, its last column (representing
    ## the response variable) must be all zeros.
    ##
    ## @item Wilkinson formula @tab A character vector of the form
    ## @qcode{'y ~ terms'} describing the response and predictor terms using
    ## Wilkinson notation.  For table input, the variable to the left of
    ## @qcode{'~'} is used as the response, overriding @var{resp_input}.
    ## @end multitable
    ##
    ## @code{@var{mdl} = LinearModel (@dots{}, @var{Name}, @var{Value},
    ## @dots{})} specifies additional options using one or more
    ## @qcode{Name-Value} pair arguments as described below.
    ##
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'Intercept'} @tab A logical scalar indicating
    ## whether to include a constant (intercept) term in the model.  Default
    ## is @qcode{true}.  Ignored when @var{modelspec} is a Wilkinson formula.
    ##
    ## @item @qcode{'Weights'} @tab A numeric vector of nonnegative
    ## observation weights, with one element per observation, used to fit a
    ## weighted least squares model.  Default is a vector of ones.
    ##
    ## @item @qcode{'Exclude'} @tab A numeric or logical vector
    ## specifying observations to exclude from the fit, given as row indices
    ## or a logical mask.  Excluded observations, together with any
    ## observation containing a missing value, are recorded in
    ## @code{ObservationInfo} but do not contribute to the fit.
    ##
    ## @item @qcode{'CategoricalVars'} @tab Specifies which predictor
    ## variables are treated as categorical, given as a vector of column
    ## indices, a logical vector, or a cell array of variable names.  Each
    ## categorical predictor with @math{L} categories is expanded into
    ## @math{L-1} indicator (dummy) variables, using the first category as
    ## the reference level.
    ##
    ## @item @qcode{'VarNames'} @tab A cell array of character vectors
    ## naming the predictor and response variables, in order, with the
    ## response variable name last.  Only applies to matrix input, since
    ## table variables already carry their own names.
    ##
    ## @item @qcode{'ResponseVar'} @tab A character vector naming the
    ## response variable, used to override the response variable name that
    ## would otherwise be used.
    ##
    ## @item @qcode{'PredictorVars'} @tab A cell array of character
    ## vectors naming which variables in @var{tbl} to use as predictors.  By
    ## default, all variables other than the response variable are used.
    ##
    ## @item @qcode{'RobustOpts'} @tab Selects ordinary least squares or
    ## robust regression fitting.  This value can be @qcode{'off'} (default,
    ## ordinary least squares), @qcode{'on'} (robust fitting using the
    ## @qcode{'bisquare'} weighting function), the name of one of the
    ## weighting functions below, a function handle for a custom weighting
    ## function, or a scalar structure with fields @qcode{RobustWgtFun} and
    ## @qcode{Tune} specifying the weighting function and its tuning
    ## constant.  Robust fitting uses Iteratively Reweighted Least Squares
    ## (IRLS), refitting the model with updated observation weights until the
    ## coefficients converge.  Supported weighting function names:
    ## @qcode{'andrews'}, @qcode{'bisquare'}, @qcode{'cauchy'},
    ## @qcode{'fair'}, @qcode{'huber'}, @qcode{'logistic'}, @qcode{'ols'},
    ## @qcode{'talwar'}, @qcode{'welsch'}, each with its own default tuning
    ## constant.
    ## @end multitable
    ##
    ## @var{mdl} is returned as a @code{LinearModel} object.  If
    ## @qcode{'RobustOpts'} is anything other than @qcode{'off'}, the returned
    ## model is a robust fit rather than an ordinary least squares fit.
    ##
    ## @end deftypefn
    function this = LinearModel (varargin)
      ##   LinearModel (X, y, modelspec, NV...)
      ##   LinearModel (tbl, resp_input, modelspec, NV...)

      if (nargin == 0)
        return;
      endif

      data       = varargin{1};
      resp_input = varargin{2};
      modelspec  = varargin{3};
      nv_args    = varargin(4:end);

      opts       = lm_parse_nv (nv_args);
      is_formula = ischar (modelspec) && any (modelspec == '~');

      if (! istable (data))
        X_raw   = double (data);
        n_total = size (X_raw, 1);
        p_raw   = size (X_raw, 2);
        y_full  = double (resp_input(:));

        if (! isempty (opts.VarNames))
          if (numel (opts.VarNames) != p_raw + 1)
            error ("LinearModel: VarNames must have %d elements.", p_raw + 1);
          endif
          pred_names_raw = opts.VarNames(1:p_raw);
          resp_name      = opts.VarNames{end};
        else
          pred_names_raw = arrayfun (@(k) sprintf ("x%d", k), 1:p_raw, ...
                                     'UniformOutput', false);
          resp_name      = 'y';
        endif
        if (! isempty (opts.ResponseVar))
          resp_name = opts.ResponseVar;
        endif
        var_names_all = [pred_names_raw, {resp_name}];
        n_vars        = p_raw + 1;

      else ## 'table'
        tbl           = data;
        col_names     = tbl.Properties.VariableNames;
        n_total       = height (tbl);
        n_vars        = width (tbl);
        var_names_all = col_names;

        if (! isempty (opts.ResponseVar))
          resp_name = opts.ResponseVar;
          if (isnumeric (resp_input) && ! isempty (resp_input))
            y_ext = double (resp_input(:));
          endif
        elseif (ischar (resp_input) && ! isempty (resp_input))
          resp_name = resp_input;
        elseif (isstring (resp_input) && ! isempty (resp_input))
          resp_name = char (resp_input);
        elseif (isnumeric (resp_input) && ! isempty (resp_input))
          resp_name = 'y';
          y_ext = double (resp_input(:));
        elseif (is_formula)
          tparts    = strsplit (modelspec, '~');
          resp_name = strtrim (tparts{1});
        else
          resp_name = col_names{end};
        endif

        if (opts.PredictorVarsGiven)
          pv = opts.PredictorVars;
          if (isnumeric (pv) || islogical (pv))
            pred_names_raw = col_names(pv);
          else
            pred_names_raw = pv;
          endif
        else
          pred_names_raw = col_names(! strcmp (col_names, resp_name));
        endif
        p_raw = numel (pred_names_raw);

        if (exist ('y_ext', 'var'))
          y_full = y_ext;
        else
          y_full = double (tbl.(resp_name)(:));
        endif
      endif

      ## categorical column flags
      cat_logical = false (1, p_raw);
      if (! isempty (opts.CategoricalVars))
        cv = opts.CategoricalVars;
        if (islogical (cv))
          n_cv = min (numel (cv), p_raw);
          cat_logical(1:n_cv) = cv(1:n_cv);
        elseif (isnumeric (cv))
          valid_cv = cv(cv > 0 & cv <= p_raw);
          cat_logical(valid_cv) = true;
        elseif (iscell (cv))
          for i = 1:numel (cv)
            cat_logical(strcmp (pred_names_raw, cv{i})) = true;
          endfor
        endif
      endif
      if (istable (data))
        for j = 1:p_raw
          col = tbl.(pred_names_raw{j});
          ## A logical or string column groups its observations just as a cell
          ## or categorical one does, and is coded the same way.
          if (iscell (col) || isa (col, 'categorical') ...
              || islogical (col) || isa (col, 'string'))
            cat_logical(j) = true;
          endif
        endfor
      endif

      ## missing and excluded masks
      if (! istable (data))
        missing_mask = any (isnan (X_raw), 2) | isnan (y_full);
      else
        missing_mask = any (ismissing (tbl), 2);
      endif

      excluded_mask = false (n_total, 1);
      if (! isempty (opts.Exclude))
        ex = opts.Exclude(:);
        if (islogical (ex))
          excluded_mask(1:numel (ex)) = ex;
        else
          excluded_mask(ex) = true;
        endif
      endif

      subset_mask = ! missing_mask & ! excluded_mask;
      n_obs       = sum (subset_mask);

      if (n_obs < 1)
        error ("LinearModel: No observations remain after removing missing/excluded rows.");
      endif

      ## weights
      if (isempty (opts.Weights))
        w_full = ones (n_total, 1);
      else
        w_full = double (opts.Weights(:));
      endif
      w_sub = w_full(subset_mask);

      if (is_formula)

        ## PATH A: Wilkinson formula string
        if (! istable (data))
          tbl_temp = array2table ([X_raw, y_full], 'VariableNames', var_names_all);
          tbl_sub  = tbl_temp(subset_mask, :);
        else
          tbl_sub = tbl(subset_mask, :);
        endif

        [X_design_sub, ~, coef_names_raw] = parseWilkinsonFormula ( ...
          modelspec, 'model_matrix', tbl_sub, pred_names_raw(cat_logical));

        coef_names    = coef_names_raw(:)';
        y_sub         = y_full(subset_mask);
        n_coef        = size (X_design_sub, 2);
        has_intercept = any (strcmp (coef_names, '(Intercept)'));

        [terms, cat_info, term_cols] = terms_from_coefnames (coef_names, ...
                              pred_names_raw, cat_logical, data, tbl_sub);
        enc_names     = term_cols(1:end-1);

      else

        ## PATH B: Keyword / numeric terms matrix
        X_num_full     = zeros (n_total, p_raw);
        cat_str_levels = cell (1, p_raw);

        for j = 1:p_raw
          if (istable (data))
            col = tbl.(pred_names_raw{j});
            if (iscell (col))

              ## Appearance order, so that the omitted reference level is the
              ## one the data shows first; see parseWilkinsonFormula for the
              ## formula path.
              [cat_str_levels{j}, ~, ic] = unique (col, 'stable');
              X_num_full(:, j) = ic;
            elseif (isa (col, 'categorical'))
              cat_str_levels{j} = categories (col);
              [~, ic] = ismember (cellstr (col), cat_str_levels{j});
              X_num_full(:, j) = ic;
            else
              X_num_full(:, j) = double (col(:));
              cat_str_levels{j} = {};
            endif
          else
            if (cat_logical(j))
              uvals = sort (unique (X_raw(isfinite (X_raw(:,j)), j)));
              cat_str_levels{j} = strtrim (cellstr (num2str (uvals(:))));
              [~, ic] = ismember (X_raw(:,j), uvals);
              X_num_full(:, j) = ic;
            else
              X_num_full(:, j) = X_raw(:, j);
              cat_str_levels{j} = {};
            endif
          endif
        endfor

        X_num_sub = X_num_full(subset_mask, :);
        y_sub     = y_full(subset_mask);

        ## Whether a categorical is given all its indicator columns depends
        ## on the intercept, which has to be settled before encoding.
        [X_enc_sub, enc_names, cat_info] = encode_categorical ( ...
          X_num_sub, cat_logical, pred_names_raw, cat_str_levels, ...
          modelspec_has_intercept (modelspec, opts.Intercept));
        p_enc  = size (X_enc_sub, 2);

        [terms, has_intercept, coef_names, emsg] = parse_modelspec ( ...
          modelspec, enc_names, p_enc, opts.Intercept);
        if (! isempty (emsg))
          error ("LinearModel: %s", emsg);
        endif
        n_coef = rows (terms);

        X_design_sub = build_design (terms, X_enc_sub);
        term_cols    = [enc_names, {''}];

      endif

      if (isempty (opts.RobustOpts))
        fit     = LinearModel.lm_fit (X_design_sub, y_sub, w_sub);
        RobustS = [];
      else
        fit     = lm_robust_fit (X_design_sub, y_sub, w_sub, ...
                                  opts.RobustOpts.WgtFun, opts.RobustOpts.Tune);
        RobustS.RobustWgtFun = opts.RobustOpts.WgtFun;
        RobustS.Tune         = opts.RobustOpts.Tune;
        RobustS.Weights      = fit.RobustWeights;
      endif

      D   = lm_diagnostics (X_design_sub, y_sub, fit, w_sub);

      p    = fit.rank_X;
      SSE  = fit.SSE;
      SSR  = fit.SSR;
      SST  = fit.SST;
      DFE  = fit.DFE;
      MSE  = fit.MSE;
      RMSE = fit.RMSE;

      crit          = LinearModel.lm_criteria (fit, n_obs, has_intercept);
      LogLikelihood = crit.LogLikelihood;
      AIC           = crit.AIC;
      AICc          = crit.AICc;
      BIC           = crit.BIC;
      CAIC          = crit.CAIC;
      R2_ord        = crit.Rsquared;
      R2_adj        = crit.AdjRsquared;
      Fstat         = crit.Fstat;
      Fpval         = crit.Fpval;

      h        = fit.leverage;
      S2_i_sub = D.S2_i;

      Fitted_full = NaN (n_total, 1);
      Fitted_full(subset_mask) = fit.Fitted;

      Raw_full = NaN (n_total, 1);
      Raw_full(subset_mask) = fit.Raw;

      Pearson_sub = fit.Raw / sqrt (max (MSE, eps));
      Std_sub     = fit.Raw ./ (RMSE .* sqrt (max (1 - h, eps)));
      Stu_sub     = fit.Raw ./ (sqrt (max (S2_i_sub, eps)) .* sqrt (max (1 - h, eps)));

      Pearson_full = NaN (n_total, 1);
      Std_full     = NaN (n_total, 1);
      Stu_full     = NaN (n_total, 1);
      Pearson_full(subset_mask) = Pearson_sub;
      Std_full(subset_mask)     = Std_sub;
      Stu_full(subset_mask)     = Stu_sub;

      beta_full  = fit.beta;
      se_full    = zeros (n_coef, 1);
      tstat_full = NaN (n_coef, 1);
      pval_full  = NaN (n_coef, 1);
      active     = fit.active_cols;
      cov_diag   = diag (fit.CovBeta);

      se_full(active)    = sqrt (cov_diag(active));
      tstat_full(active) = beta_full(active) ./ se_full(active);
      pval_full(active)  = 2 * tcdf (-abs (tstat_full(active)), DFE);

      CoeffTable = table (beta_full, se_full, tstat_full, pval_full, ...
        'VariableNames', {'Estimate', 'SE', 'tStat', 'pValue'}, ...
        'RowNames',      coef_names(:));

      ResidTable = table (Raw_full, Pearson_full, Stu_full, Std_full, ...
        'VariableNames', {'Raw', 'Pearson', 'Studentized', 'Standardized'});

      Lev_full = zeros (n_total, 1);
      CD_full  = NaN (n_total, 1);
      Dff_full = NaN (n_total, 1);
      S2i_full = NaN (n_total, 1);
      CR_full  = NaN (n_total, 1);
      Lev_full(subset_mask) = D.Leverage;
      CD_full(subset_mask)  = D.CooksDistance;
      Dff_full(subset_mask) = D.Dffits;
      S2i_full(subset_mask) = D.S2_i;
      CR_full(subset_mask)  = D.CovRatio;

      Dfb_full   = NaN (n_total, n_coef);
      Dfb_full(subset_mask, :) = D.Dfbetas;

      HatMat_pad = zeros (n_total, n_total);
      HatMat_pad(subset_mask, subset_mask) = D.HatMatrix;

      DiagTable = table (Lev_full, CD_full, Dff_full, S2i_full, CR_full, ...
        Dfb_full, HatMat_pad, ...
        'VariableNames', {'Leverage', 'CooksDistance', 'Dffits', 'S2_i', ...
                          'CovRatio', 'Dfbetas', 'HatMatrix'});

      dummy_names = {};
      dummy_bases = {};
      for ci = 1:numel (cat_info.names)
        base_nm  = cat_info.names{ci};
        levels_c = cat_info.levels{ci};
        for L = 1:numel (levels_c)
          dummy_names{end+1} = [base_nm, '_', char(levels_c{L})];
          dummy_bases{end+1} = base_nm;
        endfor
      endfor

      if (has_intercept)
        orig_idx = find (! strcmp (coef_names, '(Intercept)'));
      else
        orig_idx = 1:numel (coef_names);
      endif
      non_int = coef_names(orig_idx);

      disp_terms = {};
      grp_cols   = {};
      for t = 1:numel (non_int)
        factors_t = strsplit (non_int{t}, ':');
        for f = 1:numel (factors_t)
          idx = find (strcmp (dummy_names, factors_t{f}), 1);
          if (! isempty (idx))
            factors_t{f} = dummy_bases{idx};
          endif
        endfor
        nm = strjoin (factors_t, ':');
        k = find (strcmp (disp_terms, nm), 1);
        if (isempty (k))
          disp_terms{end+1} = nm;
          grp_cols{end+1}    = orig_idx(t);
        else
          grp_cols{k}(end+1) = orig_idx(t);
        endif
      endfor

      ## Model formula.  TERMS is expressed over the encoded design columns; the
      ## formula is expressed over the model's variables, so a categorical's
      ## indicator columns have to be folded back onto the variable they came
      ## from before the term names can be built.
      var_idx = zeros (1, p_raw);
      for j = 1:p_raw
        k = find (strcmp (var_names_all, pred_names_raw{j}), 1);
        if (! isempty (k))
          var_idx(j) = k;
        endif
      endfor

      ## Indexed by the columns of TERMS, which are not the coefficient
      ## names: a factor appearing only inside an interaction or a power has a
      ## column without ever being a coefficient.
      [enc2raw, col_pow] = encodednames_to_row (term_cols(1:end-1), ...
                                             pred_names_raw, cat_info);
      terms_var = variable_level_terms (terms(:, 1:end-1), enc2raw, var_idx, ...
                                   n_vars, col_pow);

      ## A response passed apart from the table is still a variable of the
      ## model, and the formula is written over every one of them, so it takes
      ## a column of the terms matrix like any other.
      form_vars  = var_names_all(:)';
      form_terms = terms_var;
      if (! any (strcmp (form_vars, resp_name)))
        form_vars{end+1}     = resp_name;
        form_terms(:, end+1) = 0;
      endif

      FormulaObj = LinearFormula (form_terms, form_vars, ...
                                  'ResponseName', resp_name);

      ObsInfo = table (w_full, excluded_mask, missing_mask, subset_mask, ...
        'VariableNames', {'Weights', 'Excluded', 'Missing', 'Subset'});

      if (! istable (data))
        VarsTable = array2table ([X_raw, y_full], 'VariableNames', var_names_all);
      else
        VarsTable = tbl;
      endif

      nv_total   = numel (var_names_all);
      vi_class   = cell (nv_total, 1);
      vi_range   = cell (nv_total, 1);
      vi_inmodel = false (nv_total, 1);
      vi_iscat   = false (nv_total, 1);

      for j = 1:nv_total
        vname       = var_names_all{j};
        is_resp_var = strcmp (vname, resp_name);
        j_pred      = find (strcmp (pred_names_raw, vname), 1);

        if (! istable (data))
          if (! is_resp_var && ! isempty (j_pred))
            col_d = X_raw(:, j_pred);
            vi_iscat(j) = cat_logical(j_pred);
          else
            col_d = y_full;
          endif
          vi_class{j} = 'double';
          fv = col_d(subset_mask & isfinite (col_d));
          vi_range{j} = ifelse (isempty (fv), [NaN, NaN], [min(fv), max(fv)]);
        else
          col_d = tbl.(vname);
          col_d = col_d(subset_mask);

          [vi_class{j}, vi_range{j}] = variable_class_and_range (col_d);

          if (! is_resp_var && ! isempty (j_pred))
            vi_iscat(j) = cat_logical(j_pred);
          endif
        endif

        if (! is_resp_var && ! isempty (j_pred))
          vi_inmodel(j) = true;
        endif
      endfor

      VarInfo = table (vi_class, vi_range, vi_inmodel, vi_iscat, ...
        'VariableNames', {'Class', 'Range', 'InModel', 'IsCategorical'}, ...
        'RowNames',      var_names_all(:));

      this.Coefficients             = CoeffTable;
      this.CoefficientCovariance    = fit.CovBeta;
      this.CoefficientNames         = coef_names;
      this.NumCoefficients          = n_coef;
      this.NumEstimatedCoefficients = p;
      this.DFE                      = DFE;
      this.Diagnostics              = DiagTable;
      this.Fitted                   = Fitted_full;
      this.LogLikelihood            = LogLikelihood;
      this.ModelCriterion           = struct ('AIC',  AIC, 'AICc', AICc, ...
                                             'BIC',  BIC, 'CAIC', CAIC);
      if (isnan (Fstat))
        NullModelName = NaN;
      else
        NullModelName = 'constant';
      endif
      this.ModelFitVsNullModel      = struct ('Fstat',     Fstat, ...
                                             'Pvalue',    Fpval, ...
                                             'NullModel', NullModelName);
      this.MSE                      = MSE;
      this.Residuals                = ResidTable;
      this.RMSE                     = RMSE;
      this.Rsquared                 = struct ('Ordinary', R2_ord, 'Adjusted', R2_adj);
      this.SSE                      = SSE;
      this.SSR                      = SSR;
      this.SST                      = SST;
      this.Robust                   = RobustS;
      this.Steps                    = struct ();
      this.Formula                  = FormulaObj;
      this.NumObservations          = n_obs;
      this.NumPredictors            = p_raw;
      this.NumVariables             = n_vars;
      this.ObservationInfo          = ObsInfo;
      this.ObservationNames         = {};
      this.PredictorNames           = pred_names_raw;
      this.ResponseName             = resp_name;
      this.VariableInfo             = VarInfo;
      this.VariableNames            = var_names_all;
      this.Variables                = VarsTable;

      this.DesignMatrix             = X_design_sub;
      this.ActiveCols               = fit.active_cols;
      this.HasIntercept             = has_intercept;
      this.ResponseVector           = y_full;
      this.WeightVector             = w_full;
      this.SubsetMask               = subset_mask;
      this.TermsMatrix              = terms;
      this.CatLevelInfo             = cat_info;
      this.EncPredictorNames        = enc_names;
      this.EffectContrasts          = lm_effects_contrasts (this);
      this.InteractionContrasts     = lm_interaction_contrasts (this);
      this.TermGroups               = struct ('Name', disp_terms, 'Cols', grp_cols);
      this.OrigOpts                 = opts;
      if (! is_formula)
        this.EncodedPredMatrix      = X_enc_sub;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {@var{ypred} =} predict (@var{mdl}, @var{Xnew})
    ## @deftypefnx {LinearModel} {@var{ypred} =} predict (@var{mdl})
    ## @deftypefnx {LinearModel} {[@var{ypred}, @var{yci}] =} predict (@var{mdl}, @var{Xnew})
    ## @deftypefnx {LinearModel} {[@var{ypred}, @var{yci}] =} predict (@var{mdl}, @var{Xnew}, @var{Name}, @var{Value})
    ##
    ## Predict responses from a fitted linear regression model.
    ##
    ## @code{@var{ypred} = predict (@var{mdl}, @var{Xnew})} returns the fitted
    ## response values at the new predictor locations in @var{Xnew}.  @var{Xnew}
    ## can be a numeric matrix with one column per predictor in the same order
    ## as the training data, or a table whose column names match
    ## @code{@var{mdl}.PredictorNames}.  Rows containing @code{NaN} are returned
    ## as @code{NaN} without error.
    ##
    ## @code{@var{ypred} = predict (@var{mdl})} omits @var{Xnew} and returns
    ## fitted values for the original training observations in their original
    ## row order.  Rows that were excluded or contained missing values are
    ## returned as @code{NaN}.  The result is identical to
    ## @code{@var{mdl}.Fitted}.
    ##
    ## @code{[@var{ypred}, @var{yci}] = predict (@dots{})} also returns
    ## @var{yci}, an @math{n}-by-2 matrix of confidence bounds where column 1 is
    ## the lower bound and column 2 is the upper bound.  By default these are
    ## 95% pointwise confidence intervals on the mean response.
    ##
    ## Name-Value pair arguments:
    ##
    ## @multitable @columnfractions 0.2 0.78
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'Alpha'} @tab Significance level for the confidence
    ## interval, specified as a scalar in @math{[0,1]}.  The interval has
    ## coverage @math{100(1-\alpha)\%}.  Default is @code{0.05}, giving a 95%
    ## interval.
    ##
    ## @item @qcode{'Prediction'} @tab Type of interval to compute.
    ## @code{"curve"} (default) gives a confidence interval on the mean response
    ## @math{f(x)}.  @code{"observation"} gives a wider prediction interval for
    ## a single future observation @math{y = f(x) + \varepsilon}, which accounts
    ## for both estimation uncertainty and irreducible noise; it adds
    ## @code{@var{mdl}.MSE} to the variance before computing the half-width.
    ##
    ## @item @qcode{'Simultaneous'} @tab Logical flag controlling whether
    ## the bounds are simultaneous or pointwise.  When @code{true},
    ## Scheff@'{e}'s method is used so the entire predicted curve lies within
    ## the band with @math{100(1-\alpha)\%} confidence; these bands are always
    ## wider than pointwise ones.  Default is @code{false}.
    ## @end multitable
    ##
    ## @end deftypefn
    function [ypred, yci] = predict (mdl, Xnew, varargin)

      alpha    = 0.05;
      pred_obs = false;
      simultan = false;

      i = 1;
      while (i <= numel (varargin))
        if (strcmpi (varargin{i}, 'Alpha'))
          alpha = varargin{i+1};
          if (! isscalar (alpha) || ! isnumeric (alpha) || alpha < 0 || alpha > 1)
            error ("predict: Alpha must be a scalar in [0,1].");
          endif
          i += 2;
        elseif (strcmpi (varargin{i}, 'Prediction'))
          pred_str = lower (char (varargin{i+1}));
          if (! any (strcmp (pred_str, {'curve', 'observation'})))
            error ("predict: Prediction must be 'curve' or 'observation'.");
          endif
          pred_obs = strcmp (pred_str, 'observation');
          i += 2;
        elseif (strcmpi (varargin{i}, 'Simultaneous'))
          simultan = logical (varargin{i+1});
          i += 2;
        else
          error ("predict: unknown option '%s'.", varargin{i});
        endif
      endwhile

      if (nargin < 2 || isempty (Xnew))
        Xnew = mdl.Variables;
      endif

      pred_names = mdl.PredictorNames;
      p_raw      = mdl.NumPredictors;

      if (istable (Xnew))
        n_new = height (Xnew);
        X_raw = zeros (n_new, p_raw);
        for j = 1:p_raw
          if (! ismember (pred_names{j}, Xnew.Properties.VariableNames))
            error ("predict: Xnew table is missing predictor '%s'.", pred_names{j});
          endif
          col = Xnew.(pred_names{j});
          if (iscell (col))
            cat_idx = [];
            if (! isempty (mdl.CatLevelInfo.names))
              cat_idx = find (strcmp (mdl.CatLevelInfo.names, pred_names{j}));
            endif
            if (! isempty (cat_idx))
              levels_j = mdl.CatLevelInfo.levels{cat_idx};
              codes    = zeros (n_new, 1);
              for k = 1:numel (levels_j)
                codes(strcmp (col, levels_j{k})) = k;
              endfor
              X_raw(:, j) = codes;
            endif
          else
            X_raw(:, j) = double (col);
          endif
        endfor
      else
        X_raw = double (Xnew);
        if (columns (X_raw) != p_raw)
          error ("predict: Xnew must have %d columns.", p_raw);
        endif
        n_new = rows (X_raw);
      endif

      nan_rows     = any (isnan (X_raw), 2);
      X_enc_new    = reencode_predictors (X_raw, pred_names, mdl.CatLevelInfo, mdl.EncPredictorNames);
      X_design_new = build_design (mdl.TermsMatrix, X_enc_new);

      beta            = mdl.Coefficients.Estimate;
      ypred           = X_design_new * beta;
      ypred(nan_rows) = NaN;

      if (nargout > 1)
        CovB   = mdl.CoefficientCovariance;
        var_cv = sum ((X_design_new * CovB) .* X_design_new, 2);
        if (pred_obs)
          var_ci = var_cv + mdl.MSE;
        else
          var_ci = var_cv;
        endif
        p_est = mdl.NumEstimatedCoefficients;
        if (simultan)
          mult = sqrt (p_est * finv (1 - alpha, p_est, mdl.DFE));
        else
          mult = tinv (1 - alpha / 2, mdl.DFE);
        endif
        hw              = mult * sqrt (max (var_ci, 0));
        yci             = [ypred - hw, ypred + hw];
        yci(nan_rows,:) = NaN;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {@var{ysim} =} random (@var{mdl}, @var{Xnew})
    ##
    ## Simulate responses with random noise from a fitted linear regression
    ## model.
    ##
    ## @code{@var{ysim} = random (@var{mdl}, @var{Xnew})} computes the fitted
    ## response at each row of @var{Xnew} and then adds independent Gaussian
    ## noise to each value.  The noise is drawn from @math{N(0, \sigma^2)} where
    ## @math{\sigma^2} is the estimated error variance @code{@var{mdl}.MSE}
    ## (mean squared error of the fit).  The result is a column vector of the
    ## same length as the number of rows in @var{Xnew}.
    ##
    ## @var{Xnew} is required and must be non-empty.  It can be a numeric
    ## matrix with one column per predictor in the same order as the training
    ## data, or a table whose column names match
    ## @code{@var{mdl}.PredictorNames}.  Unlike @code{predict}, there is no
    ## no-argument form; the predictor locations must always be supplied
    ## explicitly.
    ##
    ## Because the added noise is drawn freshly on every call, two calls with
    ## the same @var{Xnew} will generally produce different output.  To get
    ## reproducible results, set the random seed with @code{rand ('state', s)}
    ## before calling @code{random}.
    ##
    ## For deterministic predictions without noise, use @code{predict} or
    ## @code{feval}.  @code{predict} also provides confidence intervals on the
    ## mean response.
    ##
    ## @end deftypefn
    function ysim = random (mdl, Xnew, varargin)
      if (nargin < 2)
        error ("random: Not enough input arguments.");
      endif
      if (nargin > 2)
        error ("random: Too many input arguments.");
      endif
      if (isempty (Xnew))
        error ("random: Xnew must have %d columns.", mdl.NumPredictors);
      endif
      ypred = predict (mdl, Xnew);
      ysim  = ypred + sqrt (mdl.MSE) .* randn (numel (ypred), 1);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {@var{ypred} =} feval (@var{mdl}, @var{X})
    ## @deftypefnx {LinearModel} {@var{ypred} =} feval (@var{mdl}, @var{x1}, @var{x2}, @dots{}, @var{xp})
    ##
    ## Predict responses of a fitted linear regression model using separate
    ## predictor inputs.
    ##
    ## @code{@var{ypred} = feval (@var{mdl}, @var{X})} accepts a single
    ## numeric matrix @var{X} with one column per predictor in the same order
    ## as the training data, or a table whose column names match
    ## @code{@var{mdl}.PredictorNames}.  The output is an @math{n}-by-1 column
    ## vector.  Rows that contain @code{NaN} in any predictor column are
    ## returned as @code{NaN}.
    ##
    ## @code{@var{ypred} = feval (@var{mdl}, @var{x1}, @var{x2}, @dots{},
    ## @var{xp})} accepts exactly @code{@var{mdl}.NumPredictors} separate
    ## arguments, one per predictor variable.  All non-scalar arguments must
    ## have the same size; a scalar argument is broadcast to that size
    ## automatically.  The output shape follows the shape of the non-scalar
    ## inputs: column vector inputs give a column vector output, row vector
    ## inputs give a row vector output, and all-scalar inputs give a scalar.
    ## This form is convenient when predictor data is already stored in separate
    ## vectors rather than a combined matrix.
    ##
    ## @code{feval} gives the same numerical predictions as @code{predict} but
    ## does not support confidence intervals.  Use @code{predict} when you also
    ## need bounds on the response.  Because a @code{LinearModel} object behaves
    ## like a function through @code{feval}, it can be passed directly to
    ## routines that accept a function handle, such as @code{fminsearch} or
    ## @code{integral}.
    ##
    ## @end deftypefn
    function ypred = feval (mdl, varargin)
      p_raw   = mdl.NumPredictors;
      n_extra = nargin - 1;

      if (n_extra < 1)
        error ("feval: Not enough input arguments.");
      endif

      if (n_extra == 1)

        Xnew = varargin{1};

        if (istable (Xnew))
          for j = 1:p_raw
            if (! ismember (mdl.PredictorNames{j}, Xnew.Properties.VariableNames))
              error (strcat ("feval: X does not contain one or more predictor", " variables needed for this model."));
            endif
          endfor
        else
          if (columns (double (Xnew)) != p_raw)
            error ("feval: Predictor data matrix must have %d columns.", p_raw);
          endif
        endif

        ypred = predict (mdl, Xnew);

      elseif (n_extra == p_raw)

        for i = 1:n_extra
          if (ischar (varargin{i}) || iscategorical (varargin{i}))
            if (iscategorical (varargin{i}))
              lvl_str = char (varargin{i});
            else
              lvl_str = varargin{i};
            endif
            ci = [];
            if (! isempty (mdl.CatLevelInfo.names))
              ci = find (strcmp (mdl.CatLevelInfo.names, mdl.PredictorNames{i}));
            endif
            if (isempty (ci))
              error ("feval: predictor '%s' is not categorical.", mdl.PredictorNames{i});
            endif
            levels_i = mdl.CatLevelInfo.levels{ci};
            code     = find (strcmp (levels_i, lvl_str), 1);
            if (isempty (code))
              code = NaN;
            endif
            varargin{i} = code;
          endif
        endfor

        ref_size = [];
        for i = 1:n_extra
          if (! isscalar (varargin{i}))
            s_i = size (varargin{i});
            if (isempty (ref_size))
              ref_size = s_i;
            elseif (! isequal (s_i, ref_size))
              error ("feval: All input arguments must be the same size.");
            endif
          endif
        endfor
        if (isempty (ref_size))
          ref_size = [1, 1];
        endif

        n_pts = prod (ref_size);
        Xmat  = zeros (n_pts, p_raw);
        for i = 1:n_extra
          ai = varargin{i};
          if (isscalar (ai))
            Xmat(:, i) = ai;
          else
            Xmat(:, i) = ai(:);
          endif
        endfor

        ypred = reshape (predict (mdl, Xmat), ref_size);

      else

        error (strcat ("feval: Incorrect number of input arguments. You must provide", " either %d separate predictor variable arguments, or one", " predictor matrix with %d columns."), p_raw, p_raw);

      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {@var{ci} =} coefCI (@var{mdl})
    ## @deftypefnx {LinearModel} {@var{ci} =} coefCI (@var{mdl}, @var{alpha})
    ##
    ## Confidence intervals for the coefficient estimates of a fitted linear
    ## regression model.
    ##
    ## @code{@var{ci} = coefCI (@var{mdl})} returns 95% confidence intervals
    ## for every coefficient in @var{mdl} using a default significance level of
    ## @code{0.05}.
    ##
    ## @code{@var{ci} = coefCI (@var{mdl}, @var{alpha})} uses the significance
    ## level @var{alpha}, a scalar in @math{[0, 1]}.  The resulting intervals
    ## have coverage @math{100(1-\alpha)\%}.  Setting @var{alpha} to @code{0}
    ## produces intervals of infinite width; setting it to @code{1} collapses
    ## each interval to the corresponding point estimate.
    ##
    ## The output @var{ci} is a @math{k}-by-2 numeric matrix where
    ## @math{k = } @code{@var{mdl}.NumCoefficients}.  Row @math{j} contains
    ## the interval for the @math{j}-th coefficient, whose name is stored in
    ## @code{@var{mdl}.CoefficientNames@{j@}}.  Column 1 is the lower bound and
    ## column 2 is the upper bound.  The midpoint of each interval equals the
    ## corresponding point estimate in @code{@var{mdl}.Coefficients.Estimate}.
    ##
    ## Intervals use the Wald method:
    ## @math{b_j \pm t_{(1-\alpha/2,\,\mathrm{DFE})}\,\mathrm{SE}(b_j)},
    ## where @math{b_j} is the coefficient estimate, @math{\mathrm{SE}(b_j)} is
    ## its standard error from @code{@var{mdl}.Coefficients.SE}, and the
    ## critical value is the @math{1-\alpha/2} quantile of the
    ## @math{t}-distribution with @code{@var{mdl}.DFE} degrees of freedom.
    ## In rank-deficient models, aliased coefficients have
    ## @math{\mathrm{SE} = 0} and their row in @var{ci} is @code{[0, 0]}.
    ##
    ## @end deftypefn
    function ci = coefCI (mdl, alpha)
      if (nargin > 2)
        error ("coefCI: Too many input arguments.");
      endif
      if (nargin < 2)
        alpha = 0.05;
      endif
      if (! isscalar (alpha))
        error (strcat ("coefCI: Invalid argument at position 2.", " Value must be a scalar."));
      endif
      if (! (alpha >= 0))
        error (strcat ("coefCI: Invalid argument at position 2.", " Value must be greater than or equal to 0."));
      endif
      if (alpha > 1)
        error (strcat ("coefCI: Invalid argument at position 2.", " Value must be less than or equal to 1."));
      endif

      t  = tinv (1 - alpha / 2, mdl.DFE);
      b  = mdl.Coefficients.Estimate;
      se = mdl.Coefficients.SE;
      ci = [b - t .* se, b + t .* se];

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {@var{p} =} coefTest (@var{mdl})
    ## @deftypefnx {LinearModel} {@var{p} =} coefTest (@var{mdl}, @var{H})
    ## @deftypefnx {LinearModel} {@var{p} =} coefTest (@var{mdl}, @var{H}, @var{C})
    ## @deftypefnx {LinearModel} {[@var{p}, @var{F}] =} coefTest (@dots{})
    ## @deftypefnx {LinearModel} {[@var{p}, @var{F}, @var{r}] =} coefTest (@dots{})
    ##
    ## Linear hypothesis test on the coefficients of a fitted linear regression
    ## model.
    ##
    ## @code{coefTest} tests whether one or more linear combinations of the
    ## fitted coefficients equal specified constants.  Each linear combination
    ## is encoded as a row of the contrast matrix @var{H}, and the right-hand
    ## side is given by @var{C}.
    ##
    ## @code{@var{p} = coefTest (@var{mdl})} performs the overall model F-test:
    ## it tests the joint null hypothesis that every coefficient except the
    ## intercept is zero.  The returned p-value matches the F-statistic line
    ## printed at the bottom of the model display.
    ##
    ## @code{@var{p} = coefTest (@var{mdl}, @var{H})} tests the null hypothesis
    ## @math{H \beta = 0}, where @math{\beta} is the full coefficient vector
    ## of length @math{k = } @code{@var{mdl}.NumCoefficients}.  @var{H} must be
    ## a full-rank numeric matrix with @math{k} columns; each row specifies one
    ## linear constraint.  To test a single coefficient, use a row vector with a
    ## @code{1} in that coefficient's position and zeros elsewhere; the
    ## resulting F-statistic equals the square of the corresponding t-statistic
    ## in @code{@var{mdl}.Coefficients}.  To test a categorical predictor that
    ## expands to multiple indicator columns, include one row per indicator in
    ## @var{H}.
    ##
    ## @code{@var{p} = coefTest (@var{mdl}, @var{H}, @var{C})} tests
    ## @math{H \beta = C} instead of zero.  @var{C} must be a numeric vector
    ## with the same number of elements as rows of @var{H}; both row and column
    ## vectors are accepted.
    ##
    ## The second output @var{F} is the value of the F-statistic:
    ## @math{F = (H\hat{\beta} - C)^\prime (H V H^\prime)^{-1}
    ## (H\hat{\beta} - C) / r}, where @math{V} is
    ## @code{@var{mdl}.CoefficientCovariance} and @math{r} is the number of
    ## rows of @var{H}.  The third output @var{r} is that numerator degrees of
    ## freedom; the denominator degrees of freedom is @code{@var{mdl}.DFE}.
    ## Under the null hypothesis @math{F} follows an @math{F(r, \mathrm{DFE})}
    ## distribution and the p-value is the upper-tail probability.  When
    ## @var{H} is rank-deficient but contains no @code{NaN}, both @var{p} and
    ## @var{F} are returned as @code{NaN} without an error.
    ##
    ## @end deftypefn
    function [p, F, r] = coefTest (mdl, varargin)
      if (nargout > 3)
        error ("coefTest: Too many output arguments.");
      endif
      if (numel (varargin) > 2)
        error ("coefTest: Too many input arguments.");
      endif

      k = mdl.NumCoefficients;

      if (numel (varargin) >= 1 && ! isempty (varargin{1}))

        H = varargin{1};
        if (! isnumeric (H))
          error ("coefTest: H must be a %d-by-%d numeric matrix.", size (H, 1), k);
        endif
        if (size (H, 2) != k)
          error ("coefTest: H must be a %d-by-%d numeric matrix.", size (H, 1), k);
        endif
        if (any (any (isnan (H))))
          error (strcat ("coefTest: H is not full rank and hypotheses", " are not consistent."));
        endif
        r = size (H, 1);

        if (numel (varargin) == 2)
          C = varargin{2};
          if (! isnumeric (C))
            error ("coefTest: C must be a numeric vector.");
          endif
          C = C(:);
          if (numel (C) != r)
            error ("coefTest: H must be a %d-by-%d numeric matrix.", numel (C), k);
          endif
        else
          C = zeros (r, 1);
        endif

      else

        if (mdl.HasIntercept && k > 1)
          H = [zeros(k-1, 1), eye(k-1)];
          r = k - 1;
        else
          H = eye (k);
          r = k;
        endif
        C = zeros (r, 1);

      endif

      b    = mdl.Coefficients.Estimate;
      V    = mdl.CoefficientCovariance;
      HVH  = H * V * H';
      Hb_c = H * b - C;
      if (rcond (HVH) < eps (class (HVH)))
        F = NaN;
        p = NaN;
      else
        F = (Hb_c' * (HVH \ Hb_c)) / r;
        p = betainc (mdl.DFE / (mdl.DFE + r * F), mdl.DFE / 2, r / 2);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {@var{p} =} dwtest (@var{mdl})
    ## @deftypefnx {LinearModel} {@var{p} =} dwtest (@var{mdl}, @var{method})
    ## @deftypefnx {LinearModel} {@var{p} =} dwtest (@var{mdl}, @var{method}, @var{tail})
    ## @deftypefnx {LinearModel} {[@var{p}, @var{DW}] =} dwtest (@dots{})
    ##
    ## Durbin-Watson test for serial autocorrelation of linear regression
    ## residuals.
    ##
    ## @code{dwtest} checks whether the raw residuals of @var{mdl} are
    ## correlated with their immediate neighbours in observation order, which
    ## would violate the independence assumption of ordinary least squares.
    ## The null hypothesis is that there is no autocorrelation.  A small
    ## p-value gives evidence against this and suggests that the residuals are
    ## not independent.  This test is most meaningful when the observations
    ## have a natural ordering, such as a time series.
    ##
    ## The test is based on the Durbin-Watson statistic
    ## @math{DW = \sum_{i=1}^{n-1}(e_{i+1}-e_i)^2 / \sum_{i=1}^{n}e_i^2},
    ## where @math{e_i} are the raw residuals of the active (non-excluded)
    ## observations.  The statistic always lies in @math{[0, 4]}: values near
    ## @math{2} indicate no autocorrelation, values well below @math{2}
    ## indicate positive autocorrelation (adjacent residuals tend to have the
    ## same sign), and values well above @math{2} indicate negative
    ## autocorrelation (adjacent residuals tend to alternate in sign).
    ##
    ## @var{method} controls how the p-value is computed and defaults to
    ## @qcode{'exact'}.  @qcode{'exact'} uses the eigenvalues of the
    ## projected differencing matrix together with Imhof's numerical
    ## integration to obtain a precise p-value; this is slower but accurate
    ## for any sample size.  @qcode{'approximate'} uses a normal approximation
    ## based on the first two moments of the DW distribution under the null;
    ## this is faster and adequate for large samples but less reliable for
    ## small ones.  The argument is case-insensitive.
    ##
    ## @var{tail} selects the alternative hypothesis and defaults to
    ## @qcode{'both'}.  @qcode{'right'} tests for positive autocorrelation
    ## (@math{DW < 2}), @qcode{'left'} tests for negative autocorrelation
    ## (@math{DW > 2}), and @qcode{'both'} tests for autocorrelation in
    ## either direction.  The one-sided p-values always satisfy
    ## @math{p_{\mathrm{right}} + p_{\mathrm{left}} = 1}, and the two-sided
    ## p-value equals @math{2\min(p_{\mathrm{right}}, p_{\mathrm{left}})}.
    ##
    ## The second output @var{DW} is the value of the Durbin-Watson statistic
    ## itself; it does not depend on @var{method} or @var{tail}.
    ##
    ## @end deftypefn
    function [p, DW] = dwtest (mdl, varargin)
      if (nargout > 2)
        error ("dwtest: Too many output arguments.");
      endif
      if (numel (varargin) > 2)
        error ("dwtest: Too many input arguments.");
      endif

      method = 'exact';
      tail   = 'both';
      if (numel (varargin) >= 1)
        method = varargin{1};
      endif
      if (numel (varargin) == 2)
        tail = varargin{2};
      endif

      if (! ischar (method) || ! ismember (lower (method), {'exact', 'approximate'}))
        error ("dwtest: The METHOD argument must be 'approximate' or 'exact'.");
      endif
      method = lower (method);
      tail   = lower (tail);

      subset = logical (mdl.ObservationInfo.Subset);
      r      = mdl.Residuals.Raw(subset);
      r      = r(:);

      [p, DW] = dwtest (r, mdl.DesignMatrix, 'Method', method, 'Tail', tail);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {@var{NewMdl} =} addTerms (@var{mdl}, @var{terms})
    ##
    ## Add terms to a fitted linear regression model.
    ##
    ## @code{addTerms} returns a new @code{LinearModel} refitted on the same
    ## data and settings as @var{mdl} with the specified @var{terms} appended
    ## to the model formula.  The original model @var{mdl} is never modified;
    ## all settings including observation weights, excluded rows, and
    ## categorical variable encodings are carried over automatically.  To
    ## update a model in place, reassign the result:
    ## @code{@var{mdl} = addTerms (@var{mdl}, @var{terms})}.
    ##
    ## @var{terms} may be a character vector in Wilkinson notation.  Use
    ## @code{'x1'} for a main effect, @code{'x1:x2'} for a two-way
    ## interaction, @code{'x1*x2'} to add both main effects and their
    ## interaction in one step, @code{'x1 + x2^2'} to add several terms at
    ## once, or @code{'1'} to add an intercept to a no-intercept model.  A
    ## bare power term @code{'x1^2'} adds @code{x1} together with
    ## @code{x1^2} (and any intermediate powers), matching the Wilkinson
    ## hierarchy convention; power notation used inside an interaction, e.g.
    ## @code{'x1:x2^2'}, adds only that exact interaction term.  All
    ## variable names must match entries in @code{@var{mdl}.PredictorNames}.
    ##
    ## @var{terms} may also be a numeric matrix of size @var{t}-by-@var{v},
    ## where @var{t} is the number of terms to add and @var{v} equals
    ## @code{@var{mdl}.NumVariables}.  Entry @code{T(i,j)} is the exponent of
    ## variable @var{j} in term @var{i}.  For example, in a model with
    ## variables @code{x1}, @code{x2}, @code{y}: @code{[0 0 0]} is the
    ## intercept, @code{[0 1 0]} is @code{x2}, @code{[1 1 0]} is
    ## @code{x1:x2}, and @code{[2 0 0]} is @code{x1^2}.  The last column
    ## (response) is always zero.  A matrix with @code{@var{mdl}.NumPredictors}
    ## columns is also accepted and is automatically padded with a trailing
    ## zero column for the response.
    ##
    ## Terms that are already present in @var{mdl} are silently skipped.  If
    ## every specified term already exists, a warning is issued and @var{mdl}
    ## is returned unchanged.  For a categorical predictor, @code{addTerms}
    ## adds the full group of indicator variables for that predictor in one
    ## step rather than adding individual indicator columns.
    ##
    ## @end deftypefn
    function NewMdl = addTerms (mdl, terms)
      if (nargin < 2)
        error ("addTerms: Not enough input arguments.");
      endif
      if (nargin > 2)
        error ("addTerms: Too many input arguments.");
      endif

      nv   = mdl.NumVariables;
      pred = mdl.PredictorNames;

      if (isnumeric (terms) || islogical (terms))

        T = double (terms);
        if (isempty (T))
          error ("addTerms: Terms matrix must have %d columns.", nv);
        endif
        if (columns (T) == nv - 1)
          T = [T, zeros(rows (T), 1)];
        endif
        if (columns (T) != nv)
          error ("addTerms: Terms matrix must have %d columns.", nv);
        endif

      elseif (ischar (terms) || isstring (terms))

        terms_str   = strtrim (char (terms));
        T           = zeros (0, nv);
        plus_tokens = strsplit (terms_str, '+');

        for ti = 1:numel (plus_tokens)
          tok = strtrim (plus_tokens{ti});
          if (isempty (tok)); continue; endif

          if (! isempty (strfind (tok, '*')))
            star_parts = cellfun (@strtrim, strsplit (tok, '*'), ...
                                  'UniformOutput', false);
            n_sp       = numel (star_parts);
            colon_toks = {};
            for mask = 1:(2^n_sp - 1)
              sub = {};
              for bit = 1:n_sp
                if (bitand (mask, 2^(bit-1)))
                  sub{end+1} = star_parts{bit};
                endif
              endfor
              colon_toks{end+1} = strjoin (sub, ':');
            endfor
          else
            colon_toks = {tok};
          endif

          for ci = 1:numel (colon_toks)
            ctok = strtrim (colon_toks{ci});
            if (strcmp (ctok, '1'))
              T = [T; zeros(1, nv)];
            else
              parts = cellfun (@strtrim, strsplit (ctok, ':'), ...
                               'UniformOutput', false);
              if (numel (parts) == 1)
                part = parts{1};
                hat  = strfind (part, '^');
                if (isempty (hat))
                  vname = part;
                  exp   = 1;
                else
                  vname = strtrim (part(1:hat(1)-1));
                  exp   = str2double (strtrim (part(hat(1)+1:end)));
                endif
                idx = find (strcmp (pred, vname));
                if (isempty (idx))
                  error ("addTerms: Unrecognized variable: '%s'.", vname);
                endif
                for k = 1:exp
                  row         = zeros (1, nv);
                  row(idx(1)) = k;
                  T           = [T; row];
                endfor
              else
                row = zeros (1, nv);
                for pi = 1:numel (parts)
                  part = parts{pi};
                  hat  = strfind (part, '^');
                  if (isempty (hat))
                    vname = part;
                    exp   = 1;
                  else
                    vname = strtrim (part(1:hat(1)-1));
                    exp   = str2double (strtrim (part(hat(1)+1:end)));
                  endif
                  idx = find (strcmp (pred, vname));
                  if (isempty (idx))
                    error ("addTerms: Unrecognized variable: '%s'.", vname);
                  endif
                  row(idx(1)) = row(idx(1)) + exp;
                endfor
                T = [T; row];
              endif
            endif
          endfor

        endfor

      else
        error (strcat ("addTerms: Model update specification must be a model formula", " character vector or string scalar, or a terms matrix"));
      endif

      cat_info = mdl.CatLevelInfo;
      ename    = mdl.EncPredictorNames;
      n_pred   = nv - 1;

      ## Every candidate predictor's encoded column name(s), whether or not
      ## it is part of the model yet: a plain predictor occupies one column,
      ## a categorical one column per non-reference level, named exactly as
      ## reencode_predictors expects to find them.
      target_names = cell (n_pred, 1);
      for j = 1:n_pred
        ci = [];
        if (! isempty (cat_info) && isfield (cat_info, 'names') ...
            && ! isempty (cat_info.names))
          ci = find (strcmp (cat_info.names, pred{j}));
        endif
        if (isempty (ci))
          target_names{j} = pred(j);
        else
          levels_j  = cat_info.levels{ci};
          lvl_names = cell (1, numel (levels_j) - 1);
          for L = 2:numel (levels_j)
            lvl_names{L-1} = sprintf ("%s_%s", pred{j}, char (levels_j{L}));
          endfor
          target_names{j} = lvl_names;
        endif
      endfor
      target_enc = [target_names{:}];
      nc_full    = numel (target_enc) + 1;

      ## Re-slot the model's current encoded terms into that full space, so
      ## a predictor with no columns yet simply stays all zero.
      existing = zeros (rows (mdl.TermsMatrix), nc_full);
      existing(:, end) = mdl.TermsMatrix(:, end);
      for c = 1:numel (ename)
        col = find (strcmp (target_enc, ename{c}), 1);
        existing(:, col) = mdl.TermsMatrix(:, c);
      endfor

      ## Expand each requested raw-predictor row into that same full space.
      new_rows = zeros (0, nc_full);
      for i = 1:rows (T)
        orig_row = T(i, 1:n_pred);
        any_cat  = false;
        cat_rows = zeros (0, nc_full);
        cont_row = zeros (1, nc_full);
        col_off  = 0;
        for j = 1:n_pred
          n_cols = numel (target_names{j});
          if (orig_row(j) != 0)
            if (n_cols > 1)
              any_cat = true;
              for k = 1:n_cols
                r              = zeros (1, nc_full);
                r(col_off + k) = 1;
                cat_rows       = [cat_rows; r];
              endfor
            else
              cont_row(col_off + 1) = orig_row(j);
            endif
          endif
          col_off = col_off + n_cols;
        endfor
        if (any_cat)
          for k = 1:rows (cat_rows)
            new_rows = [new_rows; cat_rows(k,:) + cont_row];
          endfor
        else
          new_rows = [new_rows; cont_row];
        endif
      endfor

      is_new = false (rows (new_rows), 1);
      for i = 1:rows (new_rows)
        is_new(i) = ! any (all (existing == new_rows(i,:), 2));
      endfor
      new_rows = new_rows(is_new, :);

      if (isempty (new_rows))
        warning ("addTerms: There are no new terms among the terms you specified.");
        NewMdl = mdl;
        return;
      endif

      combined = [existing; new_rows];
      int_mask = all (combined(:, 1:end-1) == 0, 2);
      body     = combined(! int_mask, :);
      n_nonzero = sum (body(:, 1:end-1) != 0, 2);
      degree    = sum (body(:, 1:end-1), 2);
      tier      = zeros (rows (body), 1);
      tier(n_nonzero == 1 & degree == 1) = 1;
      tier(n_nonzero == 2)               = 2;
      tier(n_nonzero == 1 & degree > 1)  = 3;
      bitmask = zeros (rows (body), 1);
      for i = 1:rows (body)
        bitmask(i) = sum (2 .^ (find (body(i, 1:end-1)) - 1));
      endfor
      [~, order] = sortrows ([tier, bitmask]);
      combined   = [combined(int_mask, :); body(order, :)];
      NewMdl     = lm_refit (mdl, combined);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {@var{NewMdl} =} removeTerms (@var{mdl}, @var{terms})
    ##
    ## Remove terms from a fitted linear regression model.
    ##
    ## @code{removeTerms} returns a new @code{LinearModel} refitted on the same
    ## data and settings as @var{mdl}, but with the specified @var{terms}
    ## dropped from the model formula.  The original model @var{mdl} is never
    ## modified; all settings including observation weights, excluded rows, and
    ## categorical variable encodings are carried over automatically.  To
    ## update a model in place, reassign the result:
    ## @code{@var{mdl} = removeTerms (@var{mdl}, @var{terms})}.
    ##
    ## @var{terms} may be a character vector in Wilkinson notation.  Use
    ## @code{'x2'} to remove a main effect, @code{'x1:x2'} to remove an
    ## interaction, @code{'1'} to remove the intercept, or @code{'x1 + x2^2'}
    ## to remove several terms at once.  A bare power term @code{'x1^2'}
    ## removes @code{x1} together with @code{x1^2} (and any intermediate
    ## powers), matching the Wilkinson hierarchy convention; power notation
    ## used inside an interaction, e.g. @code{'x1:x2^2'}, removes only that
    ## exact interaction term.  The star operator @code{'x1*x2'} removes the
    ## main effects @code{x1} and @code{x2} together with their interaction
    ## @code{x1:x2} in a single call, following the same expansion rule as
    ## @code{addTerms}.  All variable names must match entries in
    ## @code{@var{mdl}.PredictorNames}.
    ##
    ## @var{terms} may also be a numeric matrix of size @var{t}-by-@var{v},
    ## where @var{t} is the number of terms to remove and @var{v} equals
    ## @code{@var{mdl}.NumVariables}.  Entry @code{T(i,j)} is the exponent of
    ## variable @var{j} in term @var{i}.  For example, in a model with
    ## variables @code{x1}, @code{x2}, @code{y}: @code{[0 0 0]} is the
    ## intercept, @code{[0 1 0]} is @code{x2}, @code{[1 1 0]} is
    ## @code{x1:x2}, and @code{[2 0 0]} is @code{x1^2}.  A matrix with
    ## @code{@var{mdl}.NumPredictors} columns is also accepted and is
    ## automatically padded with a trailing zero column for the response.
    ##
    ## Terms specified but absent from @var{mdl} are silently skipped.  A
    ## warning is issued and @var{mdl} is returned unchanged only when every
    ## single specified term is absent from the model.  For a categorical
    ## predictor, @code{removeTerms} removes the full group of indicator
    ## variables for that predictor in one step.
    ##
    ## @end deftypefn
    function NewMdl = removeTerms (mdl, terms)
      if (nargin < 2)
        error ("removeTerms: Not enough input arguments.");
      endif
      if (nargin > 2)
        error ("removeTerms: Too many input arguments.");
      endif

      nv   = mdl.NumVariables;
      pred = mdl.PredictorNames;

      if (isnumeric (terms) || islogical (terms))

        T = double (terms);
        if (isempty (T))
          error ("removeTerms: Terms matrix must have %d columns.", nv);
        endif
        if (columns (T) == nv - 1)
          T = [T, zeros(rows (T), 1)];
        endif
        if (columns (T) != nv)
          error ("removeTerms: Terms matrix must have %d columns.", nv);
        endif

      elseif (ischar (terms) || isstring (terms))

        terms_str   = strtrim (char (terms));
        T           = zeros (0, nv);
        plus_tokens = strsplit (terms_str, '+');

        for ti = 1:numel (plus_tokens)
          tok = strtrim (plus_tokens{ti});
          if (isempty (tok)); continue; endif

          if (! isempty (strfind (tok, '*')))
            star_parts = cellfun (@strtrim, strsplit (tok, '*'), ...
                                  'UniformOutput', false);
            n_sp       = numel (star_parts);
            colon_toks = {};
            for mask = 1:(2^n_sp - 1)
              sub = {};
              for bit = 1:n_sp
                if (bitand (mask, 2^(bit-1)))
                  sub{end+1} = star_parts{bit};
                endif
              endfor
              colon_toks{end+1} = strjoin (sub, ':');
            endfor
          else
            colon_toks = {tok};
          endif

          for ci = 1:numel (colon_toks)
            ctok = strtrim (colon_toks{ci});
            if (strcmp (ctok, '1'))
              T = [T; zeros(1, nv)];
            else
              parts = cellfun (@strtrim, strsplit (ctok, ':'), ...
                               'UniformOutput', false);
              if (numel (parts) == 1)
                part = parts{1};
                hat  = strfind (part, '^');
                if (isempty (hat))
                  vname = part;
                  exp   = 1;
                else
                  vname = strtrim (part(1:hat(1)-1));
                  exp   = str2double (strtrim (part(hat(1)+1:end)));
                endif
                idx = find (strcmp (pred, vname));
                if (isempty (idx))
                  error ("removeTerms: Unrecognized variable: '%s'.", vname);
                endif
                for k = 1:exp
                  row         = zeros (1, nv);
                  row(idx(1)) = k;
                  T           = [T; row];
                endfor
              else
                row = zeros (1, nv);
                for pi = 1:numel (parts)
                  part = parts{pi};
                  hat  = strfind (part, '^');
                  if (isempty (hat))
                    vname = part;
                    exp   = 1;
                  else
                    vname = strtrim (part(1:hat(1)-1));
                    exp   = str2double (strtrim (part(hat(1)+1:end)));
                  endif
                  idx = find (strcmp (pred, vname));
                  if (isempty (idx))
                    error ("removeTerms: Unrecognized variable: '%s'.", vname);
                  endif
                  row(idx(1)) = row(idx(1)) + exp;
                endfor
                T = [T; row];
              endif
            endif
          endfor

        endfor

      else
        error (strcat ("removeTerms: Model update specification must be a model formula", " character vector or string scalar, or a terms matrix"));
      endif

      nc = columns (mdl.TermsMatrix);
      if (nc != nv)
        cat_info    = mdl.CatLevelInfo;
        ename       = mdl.EncPredictorNames;
        n_pred      = nv - 1;
        orig_to_enc = cell (n_pred, 1);
        for j = 1:n_pred
          ci = [];
          if (! isempty (cat_info) && isfield (cat_info, 'names') ...
              && ! isempty (cat_info.names))
            ci = find (strcmp (cat_info.names, pred{j}));
          endif
          if (isempty (ci))
            orig_to_enc{j} = find (strcmp (ename, pred{j}));
          else
            levels_j = cat_info.levels{ci};
            ecols    = [];
            for L = 2:numel (levels_j)
              lvl_name = sprintf ("%s_%s", pred{j}, char (levels_j{L}));
              k        = find (strcmp (ename, lvl_name));
              if (! isempty (k))
                ecols(end+1) = k;
              endif
            endfor
            orig_to_enc{j} = ecols;
          endif
        endfor

        T_enc = zeros (0, nc);
        for i = 1:rows (T)
          orig_row = T(i, 1:n_pred);
          any_cat  = false;
          cat_rows = zeros (0, nc);
          cont_row = zeros (1, nc);
          for j = 1:n_pred
            if (orig_row(j) != 0)
              ecols = orig_to_enc{j};
              if (numel (ecols) > 1)
                any_cat = true;
                for k = 1:numel (ecols)
                  r        = zeros (1, nc);
                  r(ecols(k)) = 1;
                  cat_rows = [cat_rows; r];
                endfor
              else
                cont_row(ecols) = orig_row(j);
              endif
            endif
          endfor
          if (any_cat)
            for k = 1:rows (cat_rows)
              T_enc = [T_enc; cat_rows(k,:) + cont_row];
            endfor
          else
            T_enc = [T_enc; cont_row];
          endif
        endfor
        T = T_enc;
      endif

      existing = mdl.TermsMatrix;
      n_exist  = rows (existing);
      n_req    = rows (T);

      found = false (n_req, 1);
      for i = 1:n_req
        found(i) = any (all (existing == T(i,:), 2));
      endfor

      if (! any (found))
        warning ("removeTerms: No specified terms appear in the model.");
        NewMdl = mdl;
        return;
      endif

      keep = true (n_exist, 1);
      for i = 1:n_req
        if (found(i))
          for j = 1:n_exist
            if (keep(j) && all (existing(j,:) == T(i,:)))
              keep(j) = false;
              break;
            endif
          endfor
        endif
      endfor

      remaining = existing(keep, :);
      int_mask  = all (remaining(:, 1:end-1) == 0, 2);
      body      = remaining(! int_mask, :);
      n_nonzero = sum (body(:, 1:end-1) != 0, 2);
      degree    = sum (body(:, 1:end-1), 2);
      tier      = zeros (rows (body), 1);
      tier(n_nonzero == 1 & degree == 1) = 1;
      tier(n_nonzero == 2)               = 2;
      tier(n_nonzero == 1 & degree > 1)  = 3;
      bitmask = zeros (rows (body), 1);
      for i = 1:rows (body)
        bitmask(i) = sum (2 .^ (find (body(i, 1:end-1)) - 1));
      endfor
      [~, order] = sortrows ([tier, bitmask]);
      remaining  = [remaining(int_mask, :); body(order, :)];
      NewMdl     = lm_refit (mdl, remaining);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {} plotResiduals (@var{mdl})
    ## @deftypefnx {LinearModel} {} plotResiduals (@var{mdl}, @var{plottype})
    ## @deftypefnx {LinearModel} {} plotResiduals (@var{mdl}, @var{plottype}, @var{Name}, @var{Value})
    ## @deftypefnx {LinearModel} {} plotResiduals (@var{ax}, @dots{})
    ## @deftypefnx {LinearModel} {@var{h} =} plotResiduals (@dots{})
    ##
    ## Plot residuals of a fitted linear regression model.
    ##
    ## @code{plotResiduals (@var{mdl})} creates a probability density histogram
    ## of the raw residuals.  Bin width follows Scott's rule
    ## @math{h = 3.5 \hat\sigma n^{-1/3}} and is rounded to a visually clean
    ## value.  The bar areas sum to 1.
    ##
    ## @code{plotResiduals (@var{mdl}, @var{plottype})} creates the type of
    ## residual plot given by @var{plottype}.  For all types except
    ## @code{"histogram"} and @code{"probability"}, the full observation vector
    ## including excluded rows is passed to the plot.  Excluded or missing rows
    ## appear as @code{NaN} in the plotted data and produce visible gaps.
    ## @var{plottype} must be one of:
    ##
    ## @table @asis
    ## @item @qcode{'histogram'} (default)
    ## Probability density histogram.  Only active observations are used.
    ## Returns one @code{patch} handle.  Accepts @code{FaceColor},
    ## @code{EdgeColor}, @code{FaceAlpha}, and @code{LineWidth} Name-Value
    ## arguments.
    ##
    ## @item @qcode{'fitted'}
    ## Residuals on the y-axis against fitted values on the x-axis.  A dotted
    ## horizontal reference line marks @math{y = 0}.  Returns two line handles:
    ## @code{h(1)} is the data scatter and @code{h(2)} is the reference line.
    ##
    ## @item @qcode{'caseorder'}
    ## Residuals on the y-axis against observation row number on the x-axis,
    ## covering all rows from 1 to @code{n_total}.  A dotted horizontal
    ## reference line marks @math{y = 0}.  Returns two line handles:
    ## @code{h(1)} is the data and @code{h(2)} is the reference line.
    ##
    ## @item @qcode{'lagged'}
    ## Each residual @math{r(t)} on the y-axis against the preceding residual
    ## @math{r(t-1)} on the x-axis.  Two dotted reference lines mark
    ## @math{y = 0} and @math{x = 0}.  Returns three line handles: @code{h(1)}
    ## is the scatter, @code{h(2)} is the horizontal reference, and @code{h(3)}
    ## is the vertical reference.
    ##
    ## @item @qcode{'probability'}
    ## Normal probability plot of the sorted active residuals produced by
    ## @code{normplot}.  Returns two handles: @code{h(1)} is the data line and
    ## @code{h(2)} is the fitted reference line produced by @code{normplot}.
    ## Name-Value arguments are not applied for this plot type.
    ##
    ## @item @qcode{'observed'}
    ## Observed response values on the y-axis against fitted values on the
    ## x-axis.  A dotted @math{y = x} reference line is drawn through the
    ## origin.  Vertical segments connect each observed point down to the
    ## reference line.
    ## Returns three handles: @code{h(1)} is the scatter, @code{h(2)} is the
    ## @math{y = x} reference, and @code{h(3)} is the vertical segment line
    ## (stored as a single @code{NaN}-separated line object).
    ##
    ## @item @qcode{'symmetry'}
    ## Upper-tail distances from the median plotted against lower-tail distances
    ## from the median.  Each point @code{(x, y)} satisfies
    ## @math{x = \mathrm{med} - r_{(i)}} and
    ## @math{y = r_{(n+1-i)} - \mathrm{med}}, using the
    ## @math{\lfloor n/2 \rfloor} most extreme observations on each side.  A
    ## perfectly symmetric distribution falls on the dotted @math{y = x}
    ## reference line.  Returns two handles: @code{h(1)} is the scatter and
    ## @code{h(2)} is the reference line.
    ## @end table
    ##
    ## @code{plotResiduals (@var{ax}, @dots{})} targets the axes object @var{ax}
    ## instead of the current axes returned by @code{gca}.
    ##
    ## @code{@var{h} = plotResiduals (@dots{})} returns a vector of graphics
    ## handles.  The number of handles depends on @var{plottype} as described
    ## above.  Name-Value arguments are applied to the data handle @code{h(1)}
    ## only.  Reference lines are always drawn with the default style and are
    ## not affected by Name-Value arguments.
    ##
    ## The following Name-Value arguments are accepted.  Arguments marked
    ## @emph{histogram only} are passed directly to the @code{patch} object and
    ## have no effect on other plot types.  Arguments marked
    ## @emph{non-histogram} are applied to the scatter marker and have no
    ## effect on the histogram.
    ##
    ## @multitable @columnfractions 0.28 0.70
    ## @headitem Name @tab Description and default
    ##
    ## @item @qcode{'ResidualType'} @tab
    ## Type of residual to plot.  One of @qcode{'raw'} (default),
    ## @qcode{'pearson'}, @qcode{'standardized'}, or @qcode{'studentized'}.
    ## Case-insensitive.  Selects the corresponding column of
    ## @code{mdl.Residuals}.
    ##
    ## @item @qcode{'Color'} @tab
    ## (@emph{non-histogram}) Marker color.
    ## Default: @code{[0.1490 0.5490 0.8660]}.
    ##
    ## @item @qcode{'Marker'} @tab
    ## (@emph{non-histogram}) Marker symbol.  Any symbol accepted by
    ## @code{plot} is valid.  Default: @qcode{'x'}.
    ##
    ## @item @qcode{'MarkerSize'} @tab
    ## (@emph{non-histogram}) Marker size in points.  Default: @code{6}.
    ##
    ## @item @qcode{'MarkerEdgeColor'} @tab
    ## (@emph{non-histogram}) Marker edge color.  Default: @qcode{'auto'}.
    ##
    ## @item @qcode{'MarkerFaceColor'} @tab
    ## (@emph{non-histogram}) Marker fill color.  Default: @qcode{'none'}.
    ##
    ## @item @qcode{'LineWidth'} @tab
    ## (@emph{non-histogram}) Width of the marker edge in points.
    ## Default: @code{0.5}.
    ##
    ## @item @qcode{'FaceColor'} @tab
    ## (@emph{histogram only}) Fill color of the histogram bars.
    ## Default: @code{[0.1490 0.5490 0.8660]}.
    ##
    ## @item @qcode{'EdgeColor'} @tab
    ## (@emph{histogram only}) Edge color of the histogram bars.
    ##
    ## @item @qcode{'FaceAlpha'} @tab
    ## (@emph{histogram only}) Transparency of the histogram bars,
    ## specified as a scalar in @math{[0, 1]}.
    ## @end multitable
    ##
    ## @end deftypefn
    function h = plotResiduals (this, varargin)
      [ax, mdl, args] = lm_plot_axes (this, varargin);

      DEF_COLOR = [0.1490, 0.5490, 0.8660];
      REF_COLOR = [0.8510, 0.8510, 0.8510];

      valid_pt  = {'histogram', 'fitted', 'lagged', 'caseorder', ...
                   'probability', 'observed', 'symmetry'};
      known_nv  = {'residualtype', 'color', 'marker', 'markersize', ...
                   'markeredgecolor', 'markerfacecolor', 'linewidth', ...
                   'facecolor', 'edgecolor', 'facealpha'};

      if (! isempty (args) && ischar (args{1}) ...
          && ! any (strcmpi (args{1}, known_nv)))
        pt_str = args{1};
        args   = args(2:end);
        idx = find (strcmpi (pt_str, valid_pt));
        if (isempty (idx))
          error ("plotResiduals: Bad residuals plot type.");
        endif
        plottype = valid_pt{idx(1)};
      else
        plottype = 'histogram';
      endif

      residtype    = 'raw';
      nv_remaining = {};
      i            = 1;
      while (i <= numel (args))
        if (ischar (args{i}) && strcmpi (args{i}, 'ResidualType'))
          if (i + 1 > numel (args))
            error ("plotResiduals: ResidualType requires a value.");
          endif
          rt_val   = lower (char (args{i+1}));
          valid_rt = {'raw', 'pearson', 'standardized', 'studentized'};
          if (! any (strcmp (rt_val, valid_rt)))
            error (strcat ("plotResiduals: invalid ResidualType '%s'.", " Valid values are: 'Raw', 'Pearson',", " 'Standardized', 'Studentized'."), args{i+1});
          endif
          residtype = rt_val;
          i         = i + 2;
        else
          nv_remaining{end+1} = args{i};
          i = i + 1;
        endif
      endwhile

      if (isempty (ax))
        ax = gca ();
      endif

      switch (residtype)
        case 'raw';          rf = 'Raw';
        case 'pearson';      rf = 'Pearson';
        case 'standardized'; rf = 'Standardized';
        case 'studentized';  rf = 'Studentized';
      endswitch
      r = mdl.Residuals.(rf);

      switch (plottype)

        case 'histogram'
          r_act = r(! isnan (r));
        n_act = numel (r_act);
        s     = std (r_act);
        if (n_act <= 1 || s == 0)
          bw = 1;
          lo = floor (min (r_act)) - 0.5;
          hi = lo + 1;
        else
          bw_raw = 3.5 * s / (n_act ^ (1/3));
          mag    = 10 ^ floor (log10 (bw_raw));
          frac   = bw_raw / mag;
          if (frac < 1.5); nice = 1;
          elseif (frac < 2.5); nice = 2;
          elseif (frac < 4); nice = 3;
          elseif (frac < 7.5); nice = 5;
          else;                nice = 10;
          endif
          bw = nice * mag;
          lo = floor (min (r_act) / bw) * bw;
          hi = ceil (max (r_act) / bw) * bw;
        endif
        n_bins  = max (1, round ((hi - lo) / bw));
        centers = lo + bw/2 : bw : lo + bw * (n_bins - 0.5);
        [counts, ~] = hist (r_act, centers);
        dens  = counts / (n_act * bw);
        left  = lo + (0:n_bins-1) * bw;
        right = left + bw;
          Xp    = [left; left; right; right];
          Yp    = [zeros(1, n_bins); dens; dens; zeros(1, n_bins)];
          h     = patch (Xp, Yp, DEF_COLOR, 'FaceColor', DEF_COLOR, ...
                         nv_remaining{:}, 'Parent', ax);
          xlabel (ax, 'Residuals');
          ylabel (ax, 'Probability density');
          title (ax, 'Histogram of residuals');

        case 'fitted'
          fit   = mdl.Fitted;
          fin   = fit(! isnan (fit));
          props = lm_plot_props (nv_remaining);
          hold (ax, 'on');
          h(1)  = lm_plot_data (ax, fit, r, props);
          h(2)  = line ([min(fin), max(fin)], [0, 0], ...
                        'LineStyle', ':', 'Color', REF_COLOR, 'Parent', ax);
          hold (ax, 'off');
          xlabel (ax, 'Fitted values');
          ylabel (ax, 'Residuals');
          title (ax, 'Plot of residuals vs. fitted values');

        case 'caseorder'
          n_tot = numel (r);
          props = lm_plot_props (nv_remaining);
          hold (ax, 'on');
          h(1)  = lm_plot_data (ax, 1:n_tot, r, props);
          h(2)  = line ([1, n_tot], [0, 0], ...
                        'LineStyle', ':', 'Color', REF_COLOR, 'Parent', ax);
          hold (ax, 'off');
          xlabel (ax, 'Row number');
          ylabel (ax, 'Residuals');
          title (ax, 'Case order plot of residuals');

        case 'lagged'
          r_x    = r(1:end-1);
          r_y    = r(2:end);
          rx_fin = r_x(! isnan (r_x));
          ry_fin = r_y(! isnan (r_y));
          props  = lm_plot_props (nv_remaining);
          hold (ax, 'on');
          h(1)   = lm_plot_data (ax, r_x, r_y, props);
          h(2)   = line ([min(rx_fin), max(rx_fin)], [0, 0], ...
                         'LineStyle', ':', 'Color', REF_COLOR, 'Parent', ax);
          h(3)   = line ([0, 0], [min(ry_fin), max(ry_fin)], ...
                         'LineStyle', ':', 'Color', REF_COLOR, 'Parent', ax);
          hold (ax, 'off');
          xlabel (ax, 'Residual(t-1)');
          ylabel (ax, 'Residual(t)');
          title (ax, 'Plot of residuals vs. lagged residuals');

        case 'probability'
          r_act = r(! isnan (r));
          r_s   = sort (r_act);
          h     = normplot (ax, r_s);
          title (ax, 'Normal probability plot of residuals');
          xlabel (ax, 'Residuals');
          ylabel (ax, 'Probability');

        case 'observed'
          fit   = mdl.Fitted;
          obs   = mdl.Variables{:, mdl.ResponseName};
          fin   = fit(! isnan (fit));
          av    = [fin(:); obs(! isnan (fit))];
          xl    = [min(av), max(av)];
          n_tot = numel (fit);
          xv    = reshape ([fit(:)'; fit(:)'; NaN(1, n_tot)], 1, []);
          yv    = reshape ([fit(:)'; obs(:)'; NaN(1, n_tot)], 1, []);
          props = lm_plot_props (nv_remaining);
          hold (ax, 'on');
          h(1)  = lm_plot_data (ax, fit, obs, props);
          h(2)  = line (xl, xl, ...
                        'LineStyle', ':', 'Color', REF_COLOR, 'Parent', ax);
          h(3)  = line (xv, yv, ...
                        'LineStyle', '-', 'Color', REF_COLOR, 'Parent', ax);
          hold (ax, 'off');
          xlabel (ax, 'Fitted values');
          ylabel (ax, 'Observed response values');
          title (ax, 'Plot of observed vs. fitted values');

        case 'symmetry'
          r_act = r(! isnan (r));
          r_s   = sort (r_act);
          med   = median (r_s);
          m     = floor (numel (r_s) / 2);
          x_sym = sort (med - r_s(1:m));
          y_sym = sort (r_s(end-m+1:end) - med);
          mx    = max ([x_sym(:); y_sym(:)]);
          props = lm_plot_props (nv_remaining);
          hold (ax, 'on');
          h(1)  = lm_plot_data (ax, x_sym, y_sym, props);
          h(2)  = line ([0, mx], [0, mx], ...
                        'LineStyle', ':', 'Color', REF_COLOR, 'Parent', ax);
          hold (ax, 'off');
          xlabel (ax, 'Lower tail');
          ylabel (ax, 'Upper tail');
          title (ax, 'Symmetry plot of residuals around their median');

      endswitch

      if (nargout == 0)
        clear h;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {} plotDiagnostics (@var{mdl})
    ## @deftypefnx {LinearModel} {} plotDiagnostics (@var{mdl}, @var{plottype})
    ## @deftypefnx {LinearModel} {} plotDiagnostics (@var{mdl}, @var{plottype}, @var{Name}, @var{Value})
    ## @deftypefnx {LinearModel} {} plotDiagnostics (@var{ax}, @dots{})
    ## @deftypefnx {LinearModel} {@var{h} =} plotDiagnostics (@dots{})
    ##
    ## Plot observation diagnostics of a fitted linear regression model.
    ##
    ## @code{plotDiagnostics (@var{mdl})} creates a case-order plot of the
    ## leverage of each observation.  The x-axis is the observation row number
    ## running from 1 to the total number of rows including any excluded rows.
    ## A dotted horizontal reference line marks the recommended threshold
    ## @math{2p/n}, where @math{p} is @code{mdl.NumCoefficients} and @math{n}
    ## is @code{mdl.NumObservations}.
    ##
    ## @code{plotDiagnostics (@var{mdl}, @var{plottype})} creates the diagnostic
    ## plot specified by @var{plottype}.  For all types except @code{"contour"},
    ## the x-axis is the row number and covers all rows including excluded ones.
    ## Excluded rows produce @code{NaN} values in the diagnostic vectors, which
    ## appear as natural gaps in the plot with no special handling required.
    ## @var{plottype} must be one of:
    ##
    ## @table @asis
    ## @item @qcode{'leverage'} (default)
    ## Leverage of each observation (@code{mdl.Diagnostics.Leverage}).
    ## One dotted horizontal reference line at @math{2p/n}.
    ## Returns two handles: @code{h(1)} is the data scatter and @code{h(2)}
    ## is the reference line.
    ##
    ## @item @qcode{'cookd'}
    ## Cook's distance for each observation
    ## (@code{mdl.Diagnostics.CooksDistance}).  One dotted reference line at
    ## @math{3 \times \mathrm{mean(CooksDistance)}}, where the mean ignores
    ## @code{NaN} values.  Returns two handles: @code{h(1)} data, @code{h(2)}
    ## reference.
    ##
    ## @item @qcode{'covratio'}
    ## Delete-1 ratio of the determinant of the coefficient covariance matrix
    ## (@code{mdl.Diagnostics.CovRatio}).  Two dotted reference lines at
    ## @math{1 - 3p/n} (lower bound) and @math{1 + 3p/n} (upper bound).
    ## Both bounds are stored as a single @code{NaN}-separated line object.
    ## Returns two handles: @code{h(1)} data, @code{h(2)} combined reference.
    ##
    ## @item @qcode{'dfbetas'}
    ## Delete-1 scaled change in each coefficient estimate
    ## (@code{mdl.Diagnostics.Dfbetas}, one column per coefficient).
    ## One line object is drawn per coefficient.  Two dotted reference lines
    ## at @math{\pm 3/\sqrt{n}} are stored as a single @code{NaN}-separated
    ## line object.  Returns @math{p+1} handles: @code{h(1)} through
    ## @code{h(p)} are the per-coefficient data lines and @code{h(p+1)} is
    ## the combined reference.  Name-Value arguments are applied to all
    ## @math{p} data handles.
    ##
    ## @item @qcode{'dffits'}
    ## Delete-1 scaled change in the fitted value
    ## (@code{mdl.Diagnostics.Dffits}).  Two dotted reference lines at
    ## @math{\pm 2\sqrt{p/n}} stored as a single @code{NaN}-separated line.
    ## Returns two handles: @code{h(1)} data, @code{h(2)} combined reference.
    ##
    ## @item @qcode{'s2_i'}
    ## Delete-1 variance estimate (@code{mdl.Diagnostics.S2_i}).  One dotted
    ## reference line at @code{mdl.MSE}.  Returns two handles: @code{h(1)}
    ## data, @code{h(2)} reference.
    ##
    ## @item @qcode{'contour'}
    ## Standardized residuals on the y-axis against leverage on the x-axis,
    ## with Cook's distance contours overlaid at levels
    ## @math{[0.05, 0.10, 0.15, 0.20, 0.25]}.  The contour surface is
    ## computed on a 31-by-30 grid over the range of the active leverage and
    ## residual values.  Returns two handles: @code{h(1)} is the data scatter
    ## (a @code{line} object) and @code{h(2)} is the contour object.
    ## @end table
    ##
    ## @code{plotDiagnostics (@var{ax}, @dots{})} targets the axes object
    ## @var{ax} instead of the current axes returned by @code{gca}.
    ##
    ## @code{@var{h} = plotDiagnostics (@dots{})} returns a vector of graphics
    ## handles.  The number of handles depends on @var{plottype} as described
    ## above.  Name-Value arguments are applied to the data handle @code{h(1)},
    ## except for @code{"dfbetas"} where they are applied to all @math{p}
    ## coefficient handles.  Reference line handles are never affected by
    ## Name-Value arguments.
    ##
    ## @multitable @columnfractions 0.28 0.70
    ## @headitem Name @tab Description and default
    ##
    ## @item @qcode{'Color'} @tab
    ## Marker color for data points.  For @code{"dfbetas"} this color is
    ## applied to all @math{p} coefficient line objects.
    ## Default: @code{[0.1490 0.5490 0.8660]}.
    ##
    ## @item @qcode{'Marker'} @tab
    ## Marker symbol.  Any symbol accepted by @code{plot} is valid.
    ## Default: @qcode{'x'}.
    ##
    ## @item @qcode{'MarkerSize'} @tab
    ## Marker size in points.  Default: @code{6}.
    ##
    ## @item @qcode{'MarkerEdgeColor'} @tab
    ## Marker edge color.  Default: @qcode{'auto'}.
    ##
    ## @item @qcode{'MarkerFaceColor'} @tab
    ## Marker fill color.  Default: @qcode{'none'}.
    ##
    ## @item @qcode{'LineWidth'} @tab
    ## Width of the marker edge in points.  Default: @code{0.5}.
    ## @end multitable
    ##
    ## @end deftypefn
    function h = plotDiagnostics (this, varargin)
      [ax, mdl, args] = lm_plot_axes (this, varargin);

      REF_COLOR = [0.8510, 0.8510, 0.8660];

      valid_pt = {'leverage', 'cookd', 'covratio', 'dfbetas', ...
                  'dffits', 's2_i', 'contour'};

      if (! isempty (args) && ischar (args{1}) ...
          && ! any (strcmpi (args{1}, ...
                    {'color','marker','markersize','markeredgecolor', ...
                     'markerfacecolor','linewidth'})))
        pt_str = args{1};
        args   = args(2:end);
        idx = find (strcmpi (pt_str, valid_pt));
        if (isempty (idx))
          error ("plotDiagnostics: Bad diagnostics plot type.");
        endif
        plottype = valid_pt{idx(1)};
      else
        plottype = 'leverage';
      endif

      props = lm_plot_props (args);

      if (isempty (ax))
        ax = gca ();
      endif

      diag_t = mdl.Diagnostics;
      p      = mdl.NumCoefficients;
      n_obs  = mdl.NumObservations;
      mse    = mdl.MSE;
      n      = numel (diag_t.Leverage);

      switch (plottype)

        case 'leverage'
          lev = diag_t.Leverage;
          ref = 2 * p / n_obs;
          hold (ax, 'on');
          h(1) = lm_plot_data (ax, 1:n, lev, props);
          h(2) = line ([0, n], [ref, ref], ...
                       'LineStyle', ':', 'Color', REF_COLOR, 'Parent', ax);
          hold (ax, 'off');
          xlabel (ax, 'Row number');
          ylabel (ax, 'Leverage');
          title (ax, 'Case order plot of leverage');

        case 'cookd'
          cd_ = diag_t.CooksDistance;
          ref = 3 * mean (cd_, 'omitnan');
          hold (ax, 'on');
          h(1) = lm_plot_data (ax, 1:n, cd_, props);
          h(2) = line ([0, n], [ref, ref], ...
                       'LineStyle', ':', 'Color', REF_COLOR, 'Parent', ax);
          hold (ax, 'off');
          xlabel (ax, 'Row number');
          ylabel (ax, 'Cook''s distance');
          title (ax, 'Case order plot of Cook''s distance');

        case 'covratio'
          cv = diag_t.CovRatio;
          lo = 1 - 3*p/n_obs;
          hi = 1 + 3*p/n_obs;
          xv = [0, n, NaN, 0, n];
          yv = [lo, lo, NaN, hi, hi];
          hold (ax, 'on');
          h(1) = lm_plot_data (ax, 1:n, cv, props);
          h(2) = line (xv, yv, ...
                       'LineStyle', ':', 'Color', REF_COLOR, 'Parent', ax);
          hold (ax, 'off');
          xlabel (ax, 'Row number');
          ylabel (ax, 'Covariance ratio');
          title (ax, 'Case order plot of covariance ratio');

        case 'dfbetas'
          db  = diag_t.Dfbetas;
          thr = 3 / sqrt (n_obs);
          xv  = [0, n, NaN, 0, n];
          yv  = [-thr, -thr, NaN, thr, thr];
          hold (ax, 'on');
          for k = 1:p
            h(k) = lm_plot_data (ax, (1:n)', db(:,k), props);
          endfor
          h(p+1) = line (xv, yv, ...
                         'LineStyle', ':', 'Color', REF_COLOR, 'Parent', ax);
          hold (ax, 'off');
          xlabel (ax, 'Row number');
          ylabel (ax, 'Scaled change in coefficients');
          title (ax, 'Case order plot of scaled change in coefficients');

        case 'dffits'
          df  = diag_t.Dffits;
          thr = 2 * sqrt (p / n_obs);
          xv  = [0, n, NaN, 0, n];
          yv  = [-thr, -thr, NaN, thr, thr];
          hold (ax, 'on');
          h(1) = lm_plot_data (ax, 1:n, df, props);
          h(2) = line (xv, yv, ...
                       'LineStyle', ':', 'Color', REF_COLOR, 'Parent', ax);
          hold (ax, 'off');
          xlabel (ax, 'Row number');
          ylabel (ax, 'Scaled change in fit');
          title (ax, 'Case order plot of scaled change in fit');

        case 's2_i'
          s2 = diag_t.S2_i;
          hold (ax, 'on');
          h(1) = lm_plot_data (ax, 1:n, s2, props);
          h(2) = line ([0, n], [mse, mse], ...
                       'LineStyle', ':', 'Color', REF_COLOR, 'Parent', ax);
          hold (ax, 'off');
          xlabel (ax, 'Row number');
          ylabel (ax, 'Leave-one-out variance');
          title (ax, 'Case order plot of leave-one-out variance');

        case 'contour'
          lev    = diag_t.Leverage;
          r_raw  = mdl.Residuals.Raw;
          act    = ! isnan (lev);
          lev_a  = lev(act);
          r_a    = r_raw(act);
          x_grid = linspace (min (lev_a), max (lev_a), 31);
          y_grid = linspace (min (r_a),   max (r_a),   30);
          [Hg, Rg] = meshgrid (x_grid, y_grid);
          Hg     = min (Hg, 1 - 1e-10);
          Zg     = Rg.^2 .* Hg ./ (p .* mse .* (1 - Hg).^2);
          levels = [0.05, 0.10, 0.15, 0.20, 0.25];
          hold (ax, 'on');
          h(1) = lm_plot_data (ax, lev, r_raw, props);
          [~, h_ct] = contour (ax, x_grid, y_grid, Zg, levels);
          h(2) = h_ct;
          hold (ax, 'off');
          xlabel (ax, 'Leverage');
          ylabel (ax, 'Residual');
          title (ax, 'Cook''s distance factorization');

      endswitch

      if (nargout == 0)
        clear h;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {} plotEffects (@var{mdl})
    ## @deftypefnx {LinearModel} {} plotEffects (@var{ax}, @var{mdl})
    ## @deftypefnx {LinearModel} {@var{h} =} plotEffects (@dots{})
    ##
    ## Plot the main effects of each predictor in a fitted linear regression
    ## model.
    ##
    ## @code{plotEffects (@var{mdl})} creates a horizontal dot-and-line plot
    ## with one row per predictor.  Each dot shows the estimated main effect on
    ## the response from changing that predictor from its minimum observed value
    ## to its maximum observed value, while holding all other predictors fixed
    ## at their observed means.  A horizontal line through each dot shows the
    ## 95% confidence interval for that effect.
    ##
    ## The main effect for predictor @var{xs} is defined as
    ## @math{g(x_{s,\max}) - g(x_{s,\min})}, where the adjusted response
    ## function @math{g} evaluates the model at the specified value of
    ## @var{xs} with all other predictors set to their observed means.
    ## For numeric predictors the sign of the effect can be positive or
    ## negative depending on the direction of the relationship.
    ##
    ## @code{plotEffects (@var{ax}, @var{mdl})} creates the plot in the axes
    ## object @var{ax} instead of the current axes returned by @code{gca}.
    ##
    ## @code{@var{h} = plotEffects (@dots{})} returns a vector of
    ## @math{p+1} graphics handles where @math{p} is the number of predictors.
    ## @code{h(1)} is the line object containing the effect estimate markers
    ## (one circle per predictor, plotted as a single line object with
    ## @code{XData} of length @math{p} and @code{YData = 1:p}).
    ## @code{h(j+1)} is the confidence interval line for predictor @math{j},
    ## with @code{XData = [ci_lo, ci_hi]} and @code{YData = [j, j]}.
    ##
    ## The y-axis tick labels follow the format
    ## @qcode{'varname: min to max'}, showing the predictor name and the
    ## minimum and maximum observed values used to compute the effect.
    ##
    ## @end deftypefn
    function h = plotEffects (this, varargin)
      [ax, mdl, args] = lm_plot_axes (this, varargin);

      if (! isempty (args))
        error ("plotEffects: Wrong number of arguments.");
      endif

      p = mdl.NumPredictors;
      if (! any (any (mdl.TermsMatrix(:, 1:end-1) != 0)))
        error ("plotEffects: Model has no predictors.");
      endif

      if (isempty (ax))
        ax = gca ();
      endif

      DEF_COLOR = [0.1490, 0.5490, 0.8660];

      pred   = mdl.PredictorNames;
      V      = mdl.CoefficientCovariance;
      beta   = mdl.Coefficients.Estimate;
      t_crit = tinv (0.975, mdl.DFE);
      cinfo  = mdl.CatLevelInfo;
      C      = mdl.EffectContrasts;

      effects = zeros (1, p);
      ci_lo   = zeros (1, p);
      ci_hi   = zeros (1, p);

      for j = 1:p
        effects(j) = C(j,:) * beta;
        SE         = sqrt (max (0, C(j,:) * V * C(j,:)'));

        ci_lo(j) = effects(j) - t_crit * SE;
        ci_hi(j) = effects(j) + t_crit * SE;
      endfor

      hold (ax, 'on');
      h(1) = plot (ax, effects, 1:p, ...
                   'LineStyle', 'none', ...
                   'Marker', 'o', ...
                   'MarkerSize', 6, ...
                   'Color', DEF_COLOR);
      for j = 1:p
        h(j+1) = line ([ci_lo(j), ci_hi(j)], [j, j], ...
                        'LineStyle', '-', ...
                        'Marker', 'none', ...
                        'Color', DEF_COLOR, ...
                        'Parent', ax);
      endfor
      hold (ax, 'off');

      rn  = mdl.VariableInfo.Properties.RowNames;
      ytl = cell (p, 1);
      for j = 1:p
        ci = [];
        if (! isempty (cinfo) && isfield (cinfo, 'names') && ! isempty (cinfo.names))
          ci = find (strcmp (cinfo.names, pred{j}));
        endif

        vidx = find (strcmp (rn, pred{j}));
        rng  = mdl.VariableInfo.Range{vidx};

        if (! isempty (ci))
          levels_j = cinfo.levels{ci};
          lo_str = char (levels_j{1});
          hi_str = char (levels_j{end});
        else
          lo_str = num2str (rng(1), '%g');
          hi_str = num2str (rng(2), '%g');
        endif
        ytl{j} = [pred{j}, ': ', lo_str, ' to ', hi_str];
      endfor

      set (ax, 'YTick', 1:p, 'YTickLabel', ytl, 'YDir', 'reverse');
      ylim  (ax, [0.5, p + 0.5]);
      xlabel (ax, 'Main Effect');
      ylabel (ax, '');
      title (ax, 'Main Effects Plot');

      if (nargout == 0)
        clear h;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {} plotAdjustedResponse (@var{mdl}, @var{var})
    ## @deftypefnx {LinearModel} {} plotAdjustedResponse (@var{mdl}, @var{var}, @var{Name}, @var{Value})
    ## @deftypefnx {LinearModel} {} plotAdjustedResponse (@var{ax}, @dots{})
    ## @deftypefnx {LinearModel} {@var{h} =} plotAdjustedResponse (@dots{})
    ##
    ## Plot the adjusted response of a fitted linear regression model against
    ## a single predictor.
    ##
    ## @code{plotAdjustedResponse (@var{mdl}, @var{var})} creates an adjusted
    ## response plot for the predictor @var{var} in the linear regression
    ## model @var{mdl}.  @var{var} is a character vector or string naming a
    ## predictor in @code{mdl.PredictorNames}, or a positive integer indexing
    ## into @code{mdl.VariableNames}.
    ##
    ## An adjusted response function describes the fitted response as a
    ## function of a single predictor, with the other predictors averaged out
    ## by averaging the fitted values over the observations used in the fit.
    ## For a model @math{y_i = f (x_{1i}, x_{2i}, @dots{}, x_{pi}) + r_i}, the
    ## adjusted response function for @math{x_1} is
    ## @math{g (x_1) = (1/n) \sum_{i=1}^n f (x_1, x_{2i}, x_{3i}, @dots{},
    ## x_{pi})}, where @math{n} is the number of observations used to fit the
    ## model.  The adjusted response data value for observation @math{i} is
    ## @math{\tilde y_i = g (x_{1i}) + r_i}.
    ##
    ## For a numeric predictor, the adjusted response function is evaluated
    ## on an evenly spaced grid of 100 points spanning the minimum to the
    ## maximum observed value of @var{var}.  For a categorical predictor, the
    ## adjusted response function is evaluated at each category level.
    ##
    ## Excluded or missing observations appear as @code{NaN} in the adjusted
    ## data and produce gaps in the plotted data points.
    ##
    ## @code{plotAdjustedResponse (@var{mdl}, @var{var}, @var{Name},
    ## @var{Value})} specifies additional Name-Value arguments applied to the
    ## adjusted data points (@code{h(1)}).  The following are accepted:
    ##
    ## @multitable @columnfractions 0.28 0.70
    ## @headitem Name @tab Description and default
    ##
    ## @item @qcode{'Color'} @tab
    ## Marker color.  Default: @code{[0.1490 0.5490 0.8660]}.
    ##
    ## @item @qcode{'Marker'} @tab
    ## Marker symbol.  Default: @qcode{'x'}.
    ##
    ## @item @qcode{'MarkerSize'} @tab
    ## Marker size in points.  Default: @code{6}.
    ##
    ## @item @qcode{'MarkerEdgeColor'} @tab
    ## Marker edge color.  Default: @qcode{'auto'}.
    ##
    ## @item @qcode{'MarkerFaceColor'} @tab
    ## Marker fill color.  Default: @qcode{'none'}.
    ##
    ## @item @qcode{'LineWidth'} @tab
    ## Width of the marker edge in points.  Default: @code{0.5}.
    ## @end multitable
    ##
    ## @code{plotAdjustedResponse (@var{ax}, @dots{})} plots into the axes
    ## object @var{ax} instead of the current axes returned by @code{gca}.
    ##
    ## @code{@var{h} = plotAdjustedResponse (@dots{})} returns a 2-by-1
    ## vector of line handles.  @code{h(1)} corresponds to the adjusted
    ## response data points and @code{h(2)} corresponds to the adjusted
    ## response function.  Name-Value arguments only affect @code{h(1)}.
    ##
    ## @end deftypefn
    function h = plotAdjustedResponse (this, varargin)
      [ax, mdl, args] = lm_plot_axes (this, varargin);

      if (isempty (args))
        error ("plotAdjustedResponse: Not enough input arguments.");
      endif

      var  = args{1};
      args = args(2:end);

      vnames = mdl.VariableNames;
      pred   = mdl.PredictorNames;

      if (ischar (var) || isstring (var))
        vname = char (var);
        if (isempty (find (strcmp (vnames, vname))))
          error ("plotAdjustedResponse: '%s' is not a variable for this fit.", vname);
        endif
      elseif (isnumeric (var) && isscalar (var))
        if (var != fix (var) || var < 1)
          error (strcat ("plotAdjustedResponse: Variable must be specified as a", " name or a positive integer."));
        endif
        if (var > numel (vnames))
          error ("plotAdjustedResponse: This model only contains %d variables.", ...
                 numel (vnames));
        endif
        vname = vnames{var};
      else
        error (strcat ("plotAdjustedResponse: Variable must be specified as a", " name or a positive integer."));
      endif

      if (strcmp (vname, mdl.ResponseName))
        error ("plotAdjustedResponse: The variable '%s' is the response in this model.", ...
               vname);
      endif

      j     = find (strcmp (pred, vname));
      pname = pred{j};

      act   = mdl.ObservationInfo.Subset;
      cinfo = mdl.CatLevelInfo;
      n_act = sum (act);
      p     = numel (pred);

      [X_act, is_cat, cat_lvls] = lm_encode_active_predictors (mdl, act, pred, cinfo);

      props = lm_plot_props (args);

      if (isempty (ax))
        ax = gca ();
      endif

      FIT_COLOR = [0.9600, 0.4660, 0.1600];
      n_total   = numel (act);
      beta      = mdl.Coefficients.Estimate;
      resid     = mdl.Residuals.Raw;

      xdata = NaN (n_total, 1);
      ydata = NaN (n_total, 1);
      xdata(act) = X_act(:, j);

      if (is_cat(j))
        levels = cat_lvls{j};
        n_lvl  = numel (levels);
        fit_y  = zeros (n_lvl, 1);
        for L = 1:n_lvl
          X_rows      = X_act;
          X_rows(:,j) = L;
          X_enc = reencode_predictors (X_rows, pred, mdl.CatLevelInfo, mdl.EncPredictorNames);
          D     = build_design (mdl.TermsMatrix, X_enc);
          fit_y(L) = mean (D * beta);
        endfor
        fit_x = (1:n_lvl)';
        ydata(act) = fit_y(X_act(:, j)) + resid(act);
      else
        x_active = X_act(:, j);
        g_active = zeros (n_act, 1);
        for i = 1:n_act
          X_rows      = X_act;
          X_rows(:,j) = x_active(i);
          X_enc = reencode_predictors (X_rows, pred, mdl.CatLevelInfo, mdl.EncPredictorNames);
          D     = build_design (mdl.TermsMatrix, X_enc);
          g_active(i) = mean (D * beta);
        endfor
        ydata(act) = g_active + resid(act);

        fit_x = linspace (min (x_active), max (x_active), 100)';
        fit_y = zeros (100, 1);
        for k = 1:100
          X_rows      = X_act;
          X_rows(:,j) = fit_x(k);
          X_enc = reencode_predictors (X_rows, pred, mdl.CatLevelInfo, mdl.EncPredictorNames);
          D     = build_design (mdl.TermsMatrix, X_enc);
          fit_y(k) = mean (D * beta);
        endfor
      endif

      hold (ax, 'on');
      h(1) = lm_plot_data (ax, xdata, ydata, props);
      set (h(1), 'DisplayName', 'Adjusted data');
      h(2) = line (fit_x, fit_y, 'Color', FIT_COLOR, 'LineStyle', '-', ...
                   'Marker', 'none', 'Parent', ax, ...
                   'DisplayName', 'Adjusted fit');
      hold (ax, 'off');

      if (is_cat(j))
        set (ax, 'XTick', 1:n_lvl, 'XTickLabel', levels);
      endif

      xlabel (ax, pname);
      ylabel (ax, ['Adjusted ', mdl.ResponseName]);
      title (ax, 'Adjusted Response Plot');
      
      hleg = legend (ax, 'show');
      set (hleg, 'Location', lm_legend_corner (xdata, ydata));

      if (nargout == 0)
        clear h;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {} plotAdded (@var{mdl})
    ## @deftypefnx {LinearModel} {} plotAdded (@var{mdl}, @var{coef})
    ## @deftypefnx {LinearModel} {} plotAdded (@var{mdl}, @var{coef}, @var{Name}, @var{Value})
    ## @deftypefnx {LinearModel} {} plotAdded (@var{ax}, @dots{})
    ## @deftypefnx {LinearModel} {@var{h} =} plotAdded (@dots{})
    ##
    ## Create an added variable plot for a fitted linear regression model.
    ##
    ## @code{plotAdded (@var{mdl})} creates an added variable plot for the
    ## whole model @var{mdl} except the constant (intercept) term.
    ##
    ## @code{plotAdded (@var{mdl}, @var{coef})} creates an added variable
    ## plot for the coefficients specified by @var{coef}.  @var{coef} is a
    ## character vector or string naming a single coefficient in
    ## @code{mdl.CoefficientNames}, the name of a categorical predictor in
    ## @code{mdl.PredictorNames} (which selects that predictor's whole group
    ## of indicator coefficients), or a vector of positive integers indexing
    ## into @code{mdl.CoefficientNames}.
    ##
    ## An added variable plot, also known as a partial regression leverage
    ## plot, illustrates the incremental effect on the response of the
    ## selected terms after removing the effects of all other terms.  For a
    ## single selected predictor @math{x_1}, the response @math{y} and
    ## @math{x_1} are each fit to all other terms:
    ## @math{y_i = g_y (x_{2i}, @dots{}, x_{pi}) + r_{yi}},
    ## @math{x_{1i} = g_x (x_{2i}, @dots{}, x_{pi}) + r_{xi}}.  The adjusted
    ## values are @math{\tilde y_i = \bar y + r_{yi}} and
    ## @math{\tilde x_{1i} = \bar x_1 + r_{xi}}.  When @var{coef} selects
    ## more than one coefficient, the selected columns of the design matrix
    ## are combined into a single direction using the unit vector
    ## @math{u = \beta / \lVert \beta \rVert}, and the added variable plot is
    ## created for that combined direction.
    ##
    ## Excluded or missing observations appear as @code{NaN} in the adjusted
    ## data and produce gaps in the plotted data points.
    ##
    ## @code{plotAdded (@var{mdl}, @var{coef}, @var{Name}, @var{Value})}
    ## specifies additional Name-Value arguments applied to the adjusted data
    ## points (@code{h(1)}).  The following are accepted:
    ##
    ## @multitable @columnfractions 0.28 0.70
    ## @headitem Name @tab Description and default
    ##
    ## @item @qcode{'Color'} @tab
    ## Marker color.  Default: @code{[0.1490 0.5490 0.8660]}.
    ##
    ## @item @qcode{'Marker'} @tab
    ## Marker symbol.  Default: @qcode{'x'}.
    ##
    ## @item @qcode{'MarkerSize'} @tab
    ## Marker size in points.  Default: @code{6}.
    ##
    ## @item @qcode{'MarkerEdgeColor'} @tab
    ## Marker edge color.  Default: @qcode{'auto'}.
    ##
    ## @item @qcode{'MarkerFaceColor'} @tab
    ## Marker fill color.  Default: @qcode{'none'}.
    ##
    ## @item @qcode{'LineWidth'} @tab
    ## Width of the marker edge in points.  Default: @code{0.5}.
    ## @end multitable
    ##
    ## @code{plotAdded (@var{ax}, @dots{})} plots into the axes object
    ## @var{ax} instead of the current axes returned by @code{gca}.
    ##
    ## @code{@var{h} = plotAdded (@dots{})} returns a 3-by-1 vector of line
    ## handles.  @code{h(1)}, @code{h(2)}, and @code{h(3)} correspond to the
    ## adjusted data points, the fitted line, and the 95% confidence bounds
    ## of the fitted line, respectively.  Name-Value arguments only affect
    ## @code{h(1)}.
    ##
    ## @end deftypefn
    function h = plotAdded (this, varargin)
      [ax, mdl, args] = lm_plot_axes (this, varargin);

      cnames = mdl.CoefficientNames;
      pred   = mdl.PredictorNames;
      ncoef  = numel (cnames);

      if (isempty (args) || ! (ischar (args{1}) || isstring (args{1}) ...
          || isnumeric (args{1})))
        J = 2:ncoef;
        if (numel (J) == 1)
          label = cnames{J};
        else
          label = 'Whole Model';
        endif
      else
        coefarg = args{1};
        args    = args(2:end);

        if (ischar (coefarg) || isstring (coefarg))
          cname = char (coefarg);
          k = find (strcmp (cnames, cname));
          if (! isempty (k))
            J = k;
          else
            cinfo = mdl.CatLevelInfo;
            ci = [];
            if (! isempty (cinfo) && isfield (cinfo, 'names') ...
                && ! isempty (cinfo.names))
              ci = find (strcmp (cinfo.names, cname));
            endif
            if (isempty (ci))
              error ("plotAdded: Bad coefficient name.");
            endif
            levels_ci = cinfo.levels{ci};
            J = zeros (1, numel (levels_ci) - 1);
            for L = 2:numel (levels_ci)
              cn = sprintf ("%s_%s", cname, char (levels_ci{L}));
              J(L-1) = find (strcmp (cnames, cn));
            endfor
          endif
        else
          J = coefarg(:)';
          if (any (J != fix (J)) || any (J < 1) || any (J > ncoef))
            error ("plotAdded: Bad coefficient number.");
          endif
        endif

        if (numel (J) == 1)
          label = cnames{J};
        else
          label  = 'Specified Terms';
          cinfo  = mdl.CatLevelInfo;
          if (! isempty (cinfo) && isfield (cinfo, 'names') ...
              && ! isempty (cinfo.names))
            for ci = 1:numel (cinfo.names)
              levels_ci = cinfo.levels{ci};
              Jc = zeros (1, numel (levels_ci) - 1);
              for L = 2:numel (levels_ci)
                cn = sprintf ("%s_%s", cinfo.names{ci}, char (levels_ci{L}));
                Jc(L-1) = find (strcmp (cnames, cn));
              endfor
              if (isequal (sort (J), sort (Jc)))
                label = cinfo.names{ci};
                break;
              endif
            endfor
          endif
        endif
      endif

      if (isempty (J))
        error ("plotAdded: Bad coefficient number.");
      endif

      act   = mdl.ObservationInfo.Subset;
      cinfo = mdl.CatLevelInfo;
      n_act = sum (act);
      p     = numel (pred);

      [X_act, ~, ~] = lm_encode_active_predictors (mdl, act, pred, cinfo);

      D = reencode_predictors (X_act, pred, cinfo, mdl.EncPredictorNames);
      D = build_design (mdl.TermsMatrix, D);

      beta  = mdl.Coefficients.Estimate;
      y_act = mdl.Variables{act, mdl.ResponseName};

      Jc = setdiff (1:ncoef, J);

      if (numel (J) == 1)
        x1_act = D(:, J);
        slope  = beta(J);
      else
        u      = beta(J) / norm (beta(J));
        x1_act = D(:, J) * u;
        slope  = norm (beta(J));
      endif

      bx = D(:,Jc) \ x1_act;
      rx = x1_act - D(:,Jc) * bx;

      if (! isempty (mdl.Robust))
        w_act = mdl.Robust.Weights(act);
      else
        w_act = mdl.ObservationInfo.Weights(act);
      endif

      x1bar_w = sum (w_act .* x1_act) / sum (w_act);
      ybar_w  = sum (w_act .* y_act)  / sum (w_act);

      ry = mdl.Residuals.Raw(act) + slope * rx;

      xtilde = x1bar_w + rx;
      ytilde = ybar_w  + ry;

      intercept_avp = mean (ytilde) - slope * mean (xtilde);

      props = lm_plot_props (args);

      if (isempty (ax))
        ax = gca ();
      endif

      FIT_COLOR = [0.9600, 0.4660, 0.1600];
      n_total   = numel (act);

      xdata = NaN (n_total, 1);
      ydata = NaN (n_total, 1);
      xdata(act) = xtilde;
      ydata(act) = ytilde;

      fit_x = linspace (min (xtilde), max (xtilde), 100)';
      fit_y = intercept_avp + slope * fit_x;

      Sxx    = sum ((xtilde - mean (xtilde)) .^ 2);
      tcrit  = tinv (0.975, mdl.DFE);
      se_pred = sqrt (mdl.MSE * (1/n_act + (fit_x - mean (xtilde)).^2 / Sxx));
      halfw   = tcrit * se_pred;

      bound_x = [fit_x; NaN; fit_x];
      bound_y = [fit_y + halfw; NaN; fit_y - halfw];

      hold (ax, 'on');
      h(1) = lm_plot_data (ax, xdata, ydata, props);
      set (h(1), 'DisplayName', 'Adjusted data');
      h(2) = line (fit_x, fit_y, 'Color', FIT_COLOR, 'LineStyle', '-', ...
                   'Marker', 'none', 'Parent', ax, ...
                   'DisplayName', sprintf ("Fit: y = %g*x", slope));
      h(3) = line (bound_x, bound_y, 'Color', FIT_COLOR, 'LineStyle', ':', ...
                   'Marker', 'none', 'Parent', ax, ...
                   'DisplayName', '95% conf. bounds');
      hold (ax, 'off');

      xlabel (ax, ['Adjusted ', label]);
      ylabel (ax, ['Adjusted ', mdl.ResponseName]);
      title (ax, ['Added Variable Plot for ', label]);
      
      hleg = legend (ax, 'show');
      set (hleg, 'Location', lm_legend_corner (xdata, ydata));

      if (nargout == 0)
        clear h;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {} plot (@var{mdl})
    ## @deftypefnx {LinearModel} {} plot (@var{ax}, @var{mdl})
    ## @deftypefnx {LinearModel} {@var{h} =} plot (@dots{})
    ##
    ## Create a default diagnostic plot for a fitted linear regression model.
    ##
    ## @code{plot (@var{mdl})} creates a plot whose type depends on the
    ## number of predictors in @var{mdl}.  If @var{mdl} has two or more
    ## predictors, @code{plot} creates an added variable plot for the whole
    ## model except the constant (intercept) term, equivalent to
    ## @code{plotAdded (@var{mdl})}.  If @var{mdl} has exactly one
    ## predictor, @code{plot} creates a scatter plot of the data together
    ## with the fitted curve and its 95% confidence bounds.  If @var{mdl}
    ## has no predictors, @code{plot} creates a histogram of the residuals,
    ## equivalent to @code{plotResiduals (@var{mdl})}.
    ##
    ## For the single-predictor case, the fitted curve and confidence
    ## bounds are computed with @code{predict}, evaluated at 100 equally
    ## spaced points spanning the observed range of the predictor when the
    ## predictor is numeric, or at each level of the predictor when it is
    ## categorical.  Excluded or missing observations appear as @code{NaN}
    ## in the data and produce gaps in the plotted points.
    ##
    ## @code{plot (@var{ax}, @var{mdl})} plots into the axes object
    ## @var{ax} instead of the current axes returned by @code{gca}.
    ##
    ## @code{@var{h} = plot (@dots{})} returns a vector of graphics object
    ## handles.  For the two-or-more-predictor and no-predictor cases, see
    ## @code{plotAdded} and @code{plotResiduals}, respectively, for the
    ## meaning of @var{h}.  For the single-predictor case, @var{h}(1),
    ## @var{h}(2), and @var{h}(3) correspond to the data points, the fitted
    ## curve, and the 95% confidence bounds of the fitted curve,
    ## respectively.
    ##
    ## @end deftypefn
    function h = plot (this, varargin)
      [ax, mdl, args] = lm_plot_axes (this, varargin);

      if (! isempty (args))
        error ("plot: Too many input arguments.");
      endif

      pred      = mdl.PredictorNames;
      cinfo     = mdl.CatLevelInfo;
      enc_names = mdl.EncPredictorNames;
      enc_active = any (mdl.TermsMatrix(:, 1:end-1) != 0, 1);

      active_pred = false (1, numel (pred));
      for c = find (enc_active)
        j = find (strcmp (pred, enc_names{c}), 1);
        if (isempty (j) && ! isempty (cinfo) && isfield (cinfo, 'names'))
          for k = 1:numel (cinfo.names)
            prefix = [cinfo.names{k}, '_'];
            if (strncmp (enc_names{c}, prefix, numel (prefix)))
              j = find (strcmp (pred, cinfo.names{k}), 1);
              break;
            endif
          endfor
        endif
        if (! isempty (j))
          active_pred(j) = true;
        endif
      endfor
      n_active = sum (active_pred);

      if (n_active == 0)
        if (isempty (ax))
          h = plotResiduals (mdl);
        else
          h = plotResiduals (ax, mdl);
        endif
      elseif (n_active >= 2)
        if (isempty (ax))
          h = plotAdded (mdl);
        else
          h = plotAdded (ax, mdl);
        endif
      else
        j1    = find (active_pred, 1);
        act   = mdl.ObservationInfo.Subset;
        n_act = sum (act);

        ci = [];
        if (! isempty (cinfo) && isfield (cinfo, 'names') ...
            && ! isempty (cinfo.names))
          ci = find (strcmp (cinfo.names, pred{j1}));
        endif

        col = mdl.Variables{act, pred{j1}};
        if (! isempty (ci))
          levels_1 = cinfo.levels{ci};
          col_str = lm_col_to_str (col);
          x_act = zeros (n_act, 1);
          for L = 1:numel (levels_1)
            x_act(strcmp (col_str, char (levels_1{L}))) = L;
          endfor
          x_grid = (1:numel (levels_1))';
        else
          x_act  = double (col(:));
          x_grid = linspace (min (x_act), max (x_act), 100)';
        endif

        y_act = mdl.Variables{act, mdl.ResponseName};

        [y_fit, yci] = mdl.predict (x_grid);

        if (isempty (ax))
          ax = gca ();
        endif

        FIT_COLOR = [0.9600, 0.4660, 0.1600];
        n_total   = numel (act);

        xdata = NaN (n_total, 1);
        ydata = NaN (n_total, 1);
        xdata(act) = x_act;
        ydata(act) = y_act;

        bound_x = [x_grid; NaN; x_grid];
        bound_y = [yci(:,1); NaN; yci(:,2)];

        props = lm_plot_props ({});

        hold (ax, 'on');
        h(1) = lm_plot_data (ax, xdata, ydata, props);
        set (h(1), 'DisplayName', 'Data');
        h(2) = line (x_grid, y_fit, 'Color', FIT_COLOR, 'LineStyle', '-', ...
                     'Marker', 'none', 'Parent', ax, 'DisplayName', 'Fit');
        h(3) = line (bound_x, bound_y, 'Color', FIT_COLOR, 'LineStyle', ':', ...
                     'Marker', 'none', 'Parent', ax, ...
                     'DisplayName', '95% conf. bounds');
        hold (ax, 'off');

        if (! isempty (ci))
          set (ax, 'XTick', 1:numel (levels_1), 'XTickLabel', levels_1);
        endif

        xlabel (ax, pred{j1});
        ylabel (ax, mdl.ResponseName);
        title (ax, [mdl.ResponseName, ' vs. ', pred{j1}]);

        hleg = legend (ax, 'show');
        set (hleg, 'Location', lm_legend_corner (xdata, ydata));
      endif

      if (nargout == 0)
        clear h;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {} plotInteraction (@var{mdl}, @var{var1}, @var{var2})
    ## @deftypefnx {LinearModel} {} plotInteraction (@var{mdl}, @var{var1}, @var{var2}, @var{ptype})
    ## @deftypefnx {LinearModel} {} plotInteraction (@var{ax}, @dots{})
    ## @deftypefnx {LinearModel} {@var{h} =} plotInteraction (@dots{})
    ##
    ## Plot the interaction effects of two predictors in a fitted linear
    ## regression model.
    ##
    ## @code{plotInteraction (@var{mdl}, @var{var1}, @var{var2})} creates a
    ## plot of the main effects of @var{var1} and @var{var2} together with
    ## their conditional effects, with horizontal lines through each effect
    ## value indicating its 95% confidence interval.  @var{var1} and
    ## @var{var2} are each a character vector or string naming a variable in
    ## @code{mdl.VariableNames}, or a positive integer indexing into
    ## @code{mdl.VariableNames}; neither may name the response variable, and
    ## they must be different variables.
    ##
    ## The main effect of a predictor is the change in the adjusted response
    ## between the two predictor values that produce the minimum and maximum
    ## adjusted response, with the other predictor averaged over its own
    ## observed values row by row.  For a numeric predictor these two values
    ## are its observed minimum and maximum; for a categorical predictor
    ## every level is evaluated and the levels producing the minimum and
    ## maximum adjusted response are used, so the effect is always
    ## nonnegative.
    ##
    ## The conditional effect of @var{var1} is its effect recomputed with
    ## @var{var2} additionally held fixed at each of a small set of
    ## conditioning values, and likewise the conditional effect of
    ## @var{var2} holds @var{var1} fixed.  The conditioning values are the
    ## observed minimum, mean of the minimum and maximum, and maximum for a
    ## numeric predictor, or every level for a categorical predictor.  When
    ## the main effect and conditional effect points for a predictor do not
    ## align vertically, the model exhibits an interaction between
    ## @var{var1} and @var{var2}.
    ##
    ## @code{plotInteraction (@var{mdl}, @var{var1}, @var{var2}, @var{ptype})}
    ## selects the plot type.  @var{ptype} is @qcode{'effects'} (default), as
    ## described above, or @qcode{'predictions'}, which instead plots the
    ## adjusted response as a function of @var{var2} for each conditioning
    ## value of @var{var1} held fixed, evaluated over 101 equally spaced
    ## points spanning the observed range of @var{var2} when @var{var2} is
    ## numeric, or at each level of @var{var2} when it is categorical.
    ##
    ## @code{plotInteraction (@var{ax}, @dots{})} plots into the axes object
    ## @var{ax} instead of the current axes returned by @code{gca}.
    ##
    ## @code{@var{h} = plotInteraction (@dots{})} returns a vector of line
    ## handles.  When @var{ptype} is @qcode{'effects'}, @code{h(1)} is the
    ## marker line through the two main effect points, @code{h(2)} and
    ## @code{h(3)} are the confidence interval lines for the main effects of
    ## @var{var1} and @var{var2}, and the remaining entries are the
    ## conditional effect points and their confidence intervals, tagged
    ## @qcode{'conditional1'} for @var{var1} and @qcode{'conditional2'} for
    ## @var{var2}.  The main effect line objects are tagged @qcode{'main'}.
    ## When @var{ptype} is @qcode{'predictions'}, each entry in @var{h}
    ## corresponds to one adjusted response curve, one per conditioning
    ## value of @var{var1}.
    ##
    ## @end deftypefn
    function h = plotInteraction (this, varargin)
      [ax, mdl, args] = lm_plot_axes (this, varargin);

      if (numel (args) < 2)
        error ("plotInteraction: Not enough input arguments.");
      endif

      var1 = args{1};
      var2 = args{2};
      args = args(3:end);

      ptype = 'effects';
      if (! isempty (args) && (ischar (args{1}) || isstring (args{1})))
        ptype = lower (char (args{1}));
        args  = args(2:end);
        if (! any (strcmp (ptype, {'effects', 'predictions'})))
          error ("plotInteraction: PTYPE must be 'effects' or 'predictions'.");
        endif
      endif
      if (! isempty (args))
        error ("plotInteraction: Too many input arguments.");
      endif

      vnames = mdl.VariableNames;

      if (ischar (var1) || isstring (var1))
        v1name = char (var1);
        if (isempty (find (strcmp (vnames, v1name))))
          error ("plotInteraction: '%s' is not a variable for this fit.", v1name);
        endif
      elseif (isnumeric (var1) && isscalar (var1))
        if (var1 != fix (var1) || var1 < 1)
          error (strcat ("plotInteraction: Variable must be specified as a", " name or a positive integer."));
        endif
        if (var1 > numel (vnames))
          error ("plotInteraction: This model only contains %d variables.", numel (vnames));
        endif
        v1name = vnames{var1};
      else
        error (strcat ("plotInteraction: Variable must be specified as a", " name or a positive integer."));
      endif
      if (strcmp (v1name, mdl.ResponseName))
        error ("plotInteraction: The variable '%s' is the response in this model.", v1name);
      endif

      if (ischar (var2) || isstring (var2))
        v2name = char (var2);
        if (isempty (find (strcmp (vnames, v2name))))
          error ("plotInteraction: '%s' is not a variable for this fit.", v2name);
        endif
      elseif (isnumeric (var2) && isscalar (var2))
        if (var2 != fix (var2) || var2 < 1)
          error (strcat ("plotInteraction: Variable must be specified as a", " name or a positive integer."));
        endif
        if (var2 > numel (vnames))
          error ("plotInteraction: This model only contains %d variables.", numel (vnames));
        endif
        v2name = vnames{var2};
      else
        error (strcat ("plotInteraction: Variable must be specified as a", " name or a positive integer."));
      endif
      if (strcmp (v2name, mdl.ResponseName))
        error ("plotInteraction: The variable '%s' is the response in this model.", v2name);
      endif

      if (strcmp (v1name, v2name))
        error ("plotInteraction: VAR1 and VAR2 must be different variables.");
      endif

      pred   = mdl.PredictorNames;
      cinfo  = mdl.CatLevelInfo;
      ename  = mdl.EncPredictorNames;
      terms  = mdl.TermsMatrix;
      act    = mdl.ObservationInfo.Subset;
      n_act  = sum (act);
      p      = numel (pred);
      beta   = mdl.Coefficients.Estimate;
      V      = mdl.CoefficientCovariance;
      t_crit = tinv (0.975, mdl.DFE);

      [X_act, is_cat, cat_lvls] = lm_encode_active_predictors (mdl, act, pred, cinfo);

      j1 = find (strcmp (pred, v1name));
      j2 = find (strcmp (pred, v2name));

      if (is_cat(j1))
        levels_1 = cat_lvls{j1};
        n_lv1    = numel (levels_1);
        g_lv1    = zeros (n_lv1, 1);
        for L = 1:n_lv1
          c_row   = lm_interaction_row (X_act, j1, L, pred, cinfo, ename, terms);
          g_lv1(L) = c_row * beta;
        endfor
        [~, i_lo1] = min (g_lv1);
        [~, i_hi1] = max (g_lv1);
        lo1 = i_lo1;  hi1 = i_hi1;
        lbl1 = [v1name, ': ', char(levels_1{i_lo1}), ' to ', char(levels_1{i_hi1})];
        grid1 = (1:n_lv1)';
        grid1_lbls = cellfun (@(s) char (s), levels_1, 'UniformOutput', false);
      else
        lo1 = min (X_act(:,j1));
        hi1 = max (X_act(:,j1));
        lbl1 = [v1name, ': ', num2str(lo1), ' to ', num2str(hi1)];
        grid1 = [lo1; (lo1+hi1)/2; hi1];
        grid1_lbls = arrayfun (@(v) num2str(v,'%g'), grid1, 'UniformOutput', false);
      endif

      if (is_cat(j2))
        levels_2 = cat_lvls{j2};
        n_lv2    = numel (levels_2);
        g_lv2    = zeros (n_lv2, 1);
        for L = 1:n_lv2
          c_row   = lm_interaction_row (X_act, j2, L, pred, cinfo, ename, terms);
          g_lv2(L) = c_row * beta;
        endfor
        [~, i_lo2] = min (g_lv2);
        [~, i_hi2] = max (g_lv2);
        lo2 = i_lo2;  hi2 = i_hi2;
        lbl2 = [v2name, ': ', char(levels_2{i_lo2}), ' to ', char(levels_2{i_hi2})];
        grid2 = (1:n_lv2)';
        grid2_lbls = cellfun (@(s) char (s), levels_2, 'UniformOutput', false);
      else
        lo2 = min (X_act(:,j2));
        hi2 = max (X_act(:,j2));
        lbl2 = [v2name, ': ', num2str(lo2), ' to ', num2str(hi2)];
        grid2 = [lo2; (lo2+hi2)/2; hi2];
        grid2_lbls = arrayfun (@(v) num2str(v,'%g'), grid2, 'UniformOutput', false);
      endif

      c_hi1 = lm_interaction_row (X_act, j1, hi1, pred, cinfo, ename, terms);
      c_lo1 = lm_interaction_row (X_act, j1, lo1, pred, cinfo, ename, terms);
      eff1  = (c_hi1 - c_lo1) * beta;
      se1   = sqrt (max (0, (c_hi1 - c_lo1) * V * (c_hi1 - c_lo1)'));

      c_hi2 = lm_interaction_row (X_act, j2, hi2, pred, cinfo, ename, terms);
      c_lo2 = lm_interaction_row (X_act, j2, lo2, pred, cinfo, ename, terms);
      eff2  = (c_hi2 - c_lo2) * beta;
      se2   = sqrt (max (0, (c_hi2 - c_lo2) * V * (c_hi2 - c_lo2)'));

      n2 = numel (grid2);
      eff_c1 = zeros (n2, 1);
      se_c1  = zeros (n2, 1);
      for k = 1:n2
        c_hi = lm_interaction_row (X_act, [j1, j2], [hi1, grid2(k)], pred, cinfo, ename, terms);
        c_lo = lm_interaction_row (X_act, [j1, j2], [lo1, grid2(k)], pred, cinfo, ename, terms);
        eff_c1(k) = (c_hi - c_lo) * beta;
        se_c1(k)  = sqrt (max (0, (c_hi - c_lo) * V * (c_hi - c_lo)'));
      endfor

      n1 = numel (grid1);
      eff_c2 = zeros (n1, 1);
      se_c2  = zeros (n1, 1);
      for k = 1:n1
        c_hi = lm_interaction_row (X_act, [j2, j1], [hi2, grid1(k)], pred, cinfo, ename, terms);
        c_lo = lm_interaction_row (X_act, [j2, j1], [lo2, grid1(k)], pred, cinfo, ename, terms);
        eff_c2(k) = (c_hi - c_lo) * beta;
        se_c2(k)  = sqrt (max (0, (c_hi - c_lo) * V * (c_hi - c_lo)'));
      endfor

      if (isempty (ax))
        ax = gca ();
      endif
      cla (ax);

      DEF_COLOR = [0.1490, 0.5490, 0.8660];
      FIT_COLOR = [0.9600, 0.4660, 0.1600];

      if (strcmp (ptype, 'effects'))

        y_main1 = 1;
        y_cond1 = (2:(1+n2))';
        y_main2 = n2 + 4;
        y_cond2 = ((n2+5):(n2+4+n1))';

        hold (ax, 'on');
        line ([0, 0], [0.5, n2 + n1 + 4.5], 'LineStyle', ':', 'Marker', 'none', ...
              'Color', [0, 0, 0], 'Parent', ax);

        h(1) = plot (ax, [eff1, eff2], [y_main1, y_main2], ...
                     'LineStyle', 'none', 'Marker', 'o', 'Color', DEF_COLOR, ...
                     'Tag', 'main');
        h(2) = line ([eff1 - t_crit*se1, eff1 + t_crit*se1], [y_main1, y_main1], ...
                     'LineStyle', '-', 'Marker', 'none', 'Color', DEF_COLOR, ...
                     'Parent', ax, 'Tag', 'main');
        h(3) = line ([eff2 - t_crit*se2, eff2 + t_crit*se2], [y_main2, y_main2], ...
                     'LineStyle', '-', 'Marker', 'none', 'Color', DEF_COLOR, ...
                     'Parent', ax, 'Tag', 'main');
        h(4) = plot (ax, eff_c1, y_cond1, ...
                     'LineStyle', 'none', 'Marker', 'o', 'Color', FIT_COLOR, ...
                     'Tag', 'conditional1');
        for k = 1:n2
          h(4+k) = line ([eff_c1(k) - t_crit*se_c1(k), eff_c1(k) + t_crit*se_c1(k)], ...
                         [y_cond1(k), y_cond1(k)], ...
                         'LineStyle', '-', 'Marker', 'none', 'Color', FIT_COLOR, ...
                         'Parent', ax, 'Tag', 'conditional1');
        endfor
        h(5+n2) = plot (ax, eff_c2, y_cond2, ...
                        'LineStyle', 'none', 'Marker', 'o', 'Color', FIT_COLOR, ...
                        'Tag', 'conditional2');
        for k = 1:n1
          h(5+n2+k) = line ([eff_c2(k) - t_crit*se_c2(k), eff_c2(k) + t_crit*se_c2(k)], ...
                            [y_cond2(k), y_cond2(k)], ...
                            'LineStyle', '-', 'Marker', 'none', 'Color', FIT_COLOR, ...
                            'Parent', ax, 'Tag', 'conditional2');
        endfor
        hold (ax, 'off');

        ytl = cell (2 + n1 + n2, 1);
        ytl{1} = lbl1;
        for k = 1:n2
          ytl{1+k} = [v2name, '=', grid2_lbls{k}];
        endfor
        ytl{2+n2} = lbl2;
        for k = 1:n1
          ytl{2+n2+k} = [v1name, '=', grid1_lbls{k}];
        endfor

        set (ax, 'YTick', [y_main1; y_cond1; y_main2; y_cond2], ...
                 'YTickLabel', ytl, 'YDir', 'reverse');
        ylim  (ax, [0.5, n2 + n1 + 4.5]);
        xlabel (ax, 'Effect');
        ylabel (ax, '');
        title  (ax, ['Interaction of ', v1name, ' and ', v2name]);

      else ## 'predictions'

        if (is_cat(j2))
          x_grid2 = (1:n_lv2)';
        else
          x_grid2 = linspace (lo2, hi2, 101)';
        endif

        hold (ax, 'on');
        line (NaN, NaN, 'Color', 'none', 'Parent', ax, 'DisplayName', v1name);

        colors   = get (ax, 'ColorOrder');
        n_colors = rows (colors);

        for k = 1:n1
          y_curve = zeros (numel (x_grid2), 1);
          for m = 1:numel (x_grid2)
            c_row = lm_interaction_row (X_act, [j1, j2], [grid1(k), x_grid2(m)], ...
                                         pred, cinfo, ename, terms);
            y_curve(m) = c_row * beta;
          endfor
          h(k) = line (x_grid2, y_curve, ...
                       'Color', colors(mod(k-1, n_colors)+1, :), ...
                       'LineStyle', '-', 'Marker', 'none', 'Parent', ax, ...
                       'DisplayName', grid1_lbls{k});
        endfor
        hold (ax, 'off');

        if (is_cat(j2))
          set (ax, 'XTick', 1:n_lv2, 'XTickLabel', grid2_lbls);
        endif

        xlabel (ax, v2name);
        ylabel (ax, ['Adjusted ', mdl.ResponseName]);
        title  (ax, ['Interaction of ', v1name, ' and ', v2name]);
        legend (ax, 'show');

      endif

      if (nargout == 0)
        clear h;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {@var{cmdl} =} compact (@var{mdl})
    ##
    ## Create a compact version of a fitted linear regression model.
    ##
    ## @code{@var{cmdl} = compact (@var{mdl})} returns a
    ## @code{CompactLinearModel} object that retains the coefficient
    ## estimates, coefficient covariance, fit statistics, model formula, and
    ## fitting method information of @var{mdl}, but discards the training
    ## data and everything derived from it.  Specifically, the following
    ## properties of @var{mdl} are not carried over and are unavailable on
    ## @var{cmdl}: @code{Fitted}, @code{Residuals}, @code{Diagnostics},
    ## @code{ObservationInfo}, @code{ObservationNames}, @code{Variables},
    ## @code{Steps}, and @code{ModelFitVsNullModel}.
    ##
    ## If @var{mdl} was fit using robust regression, the @code{Robust}
    ## structure is retained on @var{cmdl} except for its @code{Weights}
    ## field, which is always emptied; @code{RobustWgtFun} and @code{Tune}
    ## are preserved unchanged.
    ##
    ## A @code{CompactLinearModel} object consumes less memory than a
    ## @code{LinearModel} object and can still be used with @code{predict},
    ## @code{feval}, @code{random}, @code{coefCI}, and @code{coefTest}, but
    ## does not support methods that require the original training data or
    ## refitting, such as @code{addTerms}, @code{removeTerms}, @code{step},
    ## and @code{dwtest}.
    ##
    ## @seealso{LinearModel, CompactLinearModel}
    ## @end deftypefn
    function CVMdl = compact (this)
      CVMdl = CompactLinearModel (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {@var{tbl} =} anova (@var{mdl})
    ## @deftypefnx {LinearModel} {@var{tbl} =} anova (@var{mdl}, @var{anovatype})
    ## @deftypefnx {LinearModel} {@var{tbl} =} anova (@var{mdl}, @qcode{"components"}, @var{sstype})
    ##
    ## Analysis of variance for a linear regression model.
    ##
    ## @code{anova (@var{mdl})} returns a table @var{tbl} with component
    ## ANOVA statistics for every term in @var{mdl} except the constant
    ## term, computed with hierarchical (@qcode{"h"}) sums of squares.  Each
    ## row gives @code{SumSq}, @code{DF}, @code{MeanSq}, @code{F}, and
    ## @code{pValue} for the corresponding term; the trailing @qcode{Error}
    ## row gives @code{SumSq = @var{mdl}.SSE}, @code{DF = @var{mdl}.DFE},
    ## @code{MeanSq = @var{mdl}.MSE}, and @code{NaN} for @code{F} and
    ## @code{pValue}.
    ##
    ## @code{anova (@var{mdl}, @var{anovatype})} selects
    ## @qcode{"components"} (default) or @qcode{"summary"}.  For
    ## @qcode{"summary"}, @var{tbl} always contains rows @qcode{Total},
    ## @qcode{Model}, and @qcode{Residual}, and additionally @qcode{. Linear}
    ## and @qcode{. Nonlinear} whenever @var{mdl} contains an interaction
    ## term or a continuous term of degree greater than 1.  @qcode{Total}
    ## reports @code{@var{mdl}.SST} with @code{DF = NumObservations - 1};
    ## @qcode{Model} reports @code{@var{mdl}.SSR} with @code{DF =
    ## NumCoefficients - HasIntercept}; @qcode{Residual} reports
    ## @code{@var{mdl}.SSE} with @code{DF = @var{mdl}.DFE}.  Whenever the
    ## data contains two or more observations sharing identical predictor
    ## values, @var{tbl} additionally contains @qcode{. Lack of fit} and
    ## @qcode{. Pure error}, splitting @qcode{Residual} into the part
    ## explained by replicated observations and the remainder.
    ##
    ## @code{anova (@var{mdl}, @qcode{"components"}, @var{sstype})} selects
    ## the sum of squares used for the component table: @code{1} (sequential,
    ## reduction from adding each term in formula order), @code{2} (reduction
    ## from adding the term to a model containing every term that does not
    ## contain it), @qcode{"h"} (default; as Type 2, but a higher-degree
    ## term in the same continuous variable, such as a squared term, is also
    ## treated as containing the lower-degree term), or @code{3} (reduction
    ## from adding the term to a model containing every other term, with
    ## categorical predictors recoded using sum-to-zero deviation contrasts
    ## instead of @var{mdl}'s reference-level coding).  Because Type 3 uses a
    ## different coding, its @qcode{Error} row can differ from @var{mdl}.SSE
    ## and @var{mdl}.DFE when @var{mdl} is missing a lower-order relative of
    ## one of its terms (e.g. an interaction without one of its main
    ## effects, or a categorical predictor fit without an intercept).
    ##
    ## @end deftypefn
    function tbl = anova (mdl, varargin)
      if (numel (varargin) > 2)
        error ("anova: too many input arguments.");
      endif

      anovatype = "components";
      if (numel (varargin) >= 1)
        anovatype = varargin{1};
        if (! ischar (anovatype) ...
            || ! any (strcmpi (anovatype, {"summary", "components"})))
          error ("anova: ANOVATYPE must be 'summary' or 'components'.");
        endif
      endif

      sstype = "h";
      if (numel (varargin) == 2)
        if (! strcmpi (anovatype, "components"))
          error ("anova: SSTYPE can only be specified with ANOVATYPE 'components'.");
        endif
        sstype  = varargin{2};
        valid_n = isnumeric (sstype) && isscalar (sstype) && any (sstype == [1, 2, 3]);
        valid_h = ischar (sstype) && strcmpi (sstype, "h");
        if (! valid_n && ! valid_h)
          error ("anova: SSTYPE must be 1, 2, or 3.");
        endif
      endif

      groups = mdl.TermGroups;
      nterm  = numel (groups);
      term_cols = {groups.Cols};
      term_name = {groups.Name};

      if (mdl.HasIntercept)
        icol = find (strcmp (mdl.CoefficientNames, '(Intercept)'));
      else
        icol = [];
      endif

      X = mdl.DesignMatrix;
      y = mdl.ResponseVector (mdl.SubsetMask);
      w = mdl.WeightVector (mdl.SubsetMask);
      all_cols = 1:columns (X);

      if (strcmpi (anovatype, "summary"))

        is_nonlinear = false (nterm, 1);
        for k = 1:nterm
          parts = strsplit (term_name{k}, ':');
          is_nonlinear(k) = (numel (parts) > 1) || ! isempty (strfind (parts{1}, '^'));
        endfor

        SumSq = [mdl.SST; mdl.SSR];
        DF    = [mdl.NumObservations - 1; mdl.NumCoefficients - mdl.HasIntercept];
        RowNm = {"Total", "Model"};

        if (any (is_nonlinear))
          lin_cols = union (icol, cell2mat (term_cols(! is_nonlinear)));
          SS_nl    = anova_delta_sse (X, y, w, lin_cols, all_cols);
          DF_nl    = numel (all_cols) - numel (lin_cols);
          SumSq    = [SumSq; SumSq(2) - SS_nl; SS_nl];
          DF       = [DF; DF(2) - DF_nl; DF_nl];
          RowNm    = [RowNm, {". Linear", ". Nonlinear"}];
        endif

        SumSq = [SumSq; mdl.SSE];
        DF    = [DF; mdl.DFE];
        RowNm = [RowNm, {"Residual"}];

        MeanSq = SumSq ./ DF;
        F      = NaN (numel (SumSq), 1);
        pValue = NaN (numel (SumSq), 1);
        for r = 2:(numel (SumSq) - 1)
          F(r)      = MeanSq(r) / mdl.MSE;
          pValue(r) = f_pvalue (mdl.DFE, DF(r), F(r));
        endfor

        pred_names = mdl.PredictorNames;
        cat_info   = mdl.CatLevelInfo;
        p_raw      = mdl.NumPredictors;
        tbl_sub    = mdl.Variables (mdl.SubsetMask, :);
        X_raw      = anova_decode_raw (tbl_sub, pred_names, cat_info);

        [~, ~, gidx] = unique (X_raw, 'rows');
        pe_ss = 0;
        pe_df = 0;
        for g = 1:max (gidx)
          rows_g = (gidx == g);
          w_g    = w(rows_g);
          y_g    = y(rows_g);
          wmean  = sum (w_g .* y_g) / sum (w_g);
          pe_ss  = pe_ss + sum (w_g .* (y_g - wmean).^2);
          pe_df  = pe_df + sum (rows_g) - 1;
        endfor

        if (pe_df > 0)
          lof_ss = mdl.SSE - pe_ss;
          lof_df = mdl.DFE - pe_df;
          lof_ms = lof_ss / lof_df;
          pe_ms  = pe_ss / pe_df;
          lof_F  = lof_ms / pe_ms;
          lof_p  = f_pvalue (pe_df, lof_df, lof_F);

          SumSq  = [SumSq; lof_ss; pe_ss];
          DF     = [DF; lof_df; pe_df];
          MeanSq = [MeanSq; lof_ms; pe_ms];
          F      = [F; lof_F; NaN];
          pValue = [pValue; lof_p; NaN];
          RowNm  = [RowNm, {". Lack of fit", ". Pure error"}];
        endif

      else

        use_seq   = isnumeric (sstype) && sstype == 1;
        use_type3 = isnumeric (sstype) && sstype == 3;
        extended  = ischar (sstype) && strcmpi (sstype, "h");

        if (use_type3)

          pred_names = mdl.PredictorNames;
          cat_info   = mdl.CatLevelInfo;
          enc_names  = mdl.EncPredictorNames;
          p_raw      = mdl.NumPredictors;

          if (! isempty (mdl.EncodedPredMatrix))
            X_enc = mdl.EncodedPredMatrix;
          else
            tbl_sub = mdl.Variables (mdl.SubsetMask, :);
            X_raw   = anova_decode_raw (tbl_sub, pred_names, cat_info);
            X_enc   = reencode_predictors (X_raw, pred_names, cat_info, enc_names);
          endif

          ## deviation (effects) coding: reference level becomes -1
          enc_pos = 0;
          for j = 1:p_raw
            ci = find (strcmp (cat_info.names, pred_names{j}));
            if (isempty (ci))
              enc_pos = enc_pos + 1;
            else
              n_lev  = numel (cat_info.levels{ci});
              block  = enc_pos + (1:(n_lev - 1));
              is_ref = 1 - sum (X_enc(:, block), 2);
              X_enc(:, block) = X_enc(:, block) - is_ref;
              enc_pos = enc_pos + n_lev - 1;
            endif
          endfor

          D_eff   = build_design (mdl.TermsMatrix, X_enc);
          Mm      = X \ D_eff;
          is_hier = (max (max (abs (D_eff - X * Mm))) < 1e-8 * max (max (abs (D_eff)))) ...
                    && (rcond (Mm) > eps);

          if (is_hier)

            Hinv = inv (Mm);
            b    = mdl.Coefficients.Estimate;
            V    = mdl.CoefficientCovariance;

            SumSq  = zeros (nterm, 1);
            DF     = zeros (nterm, 1);
            F      = zeros (nterm, 1);
            pValue = zeros (nterm, 1);
            for k = 1:nterm
              Hk        = Hinv(term_cols{k}, :);
              Hb        = Hk * b;
              HVH       = Hk * V * Hk';
              DF(k)     = numel (term_cols{k});
              F(k)      = (Hb' * (HVH \ Hb)) / DF(k);
              SumSq(k)  = F(k) * DF(k) * mdl.MSE;
              pValue(k) = f_pvalue (mdl.DFE, DF(k), F(k));
            endfor
            MeanSq = SumSq ./ DF;

            errSumSq  = mdl.SSE;
            errDF     = mdl.DFE;
            errMeanSq = mdl.MSE;

          else

            fit_eff  = LinearModel.lm_fit (D_eff, y, w, false);
            all_cols = 1:columns (D_eff);

            SumSq = zeros (nterm, 1);
            DF    = zeros (nterm, 1);
            for k = 1:nterm
              cmp_cols = setdiff (all_cols, term_cols{k});
              [SumSq(k), ~, fit_r] = anova_delta_sse ( ...
                D_eff, y, w, cmp_cols, all_cols, fit_eff);
              DF(k) = fit_eff.rank_X - fit_r.rank_X;
            endfor

            errSumSq  = fit_eff.SSE;
            errDF     = fit_eff.DFE;
            errMeanSq = fit_eff.MSE;
            MeanSq    = SumSq ./ DF;
            F         = MeanSq / errMeanSq;
            pValue    = NaN (nterm, 1);
            ok        = DF > 0;
            pValue(ok) = f_pvalue (errDF, DF(ok), F(ok));

          endif

          SumSq  = [SumSq; errSumSq];
          DF     = [DF; errDF];
          MeanSq = [MeanSq; errMeanSq];

        else

          if (! use_seq)
            contain_mx = anova_containment (groups, extended);
          endif

          SumSq = zeros (nterm, 1);
          DF    = zeros (nterm, 1);
          for k = 1:nterm
            if (use_seq)
              cmp_cols = union (icol, cell2mat (term_cols(1:k-1)));
            else
              keep = true (1, nterm);
              keep(k) = false;
              for j = 1:nterm
                if (j != k && contain_mx(k, j))
                  keep(j) = false;
                endif
              endfor
              cmp_cols = union (icol, cell2mat (term_cols(keep)));
            endif
            full_cols = union (cmp_cols, term_cols{k});
            SumSq(k)  = anova_delta_sse (X, y, w, cmp_cols, full_cols);
            DF(k)     = numel (term_cols{k});
          endfor

          MeanSq = SumSq ./ DF;
          F      = MeanSq / mdl.MSE;
          pValue = f_pvalue (mdl.DFE, DF, F);

          SumSq  = [SumSq; mdl.SSE];
          DF     = [DF; mdl.DFE];
          MeanSq = [MeanSq; mdl.MSE];

        endif

        F      = [F; NaN];
        pValue = [pValue; NaN];
        RowNm  = [term_name, {"Error"}];

      endif

      tbl = table (SumSq, DF, MeanSq, F, pValue, ...
        'VariableNames', {'SumSq', 'DF', 'MeanSq', 'F', 'pValue'}, ...
        'RowNames', RowNm(:));

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearModel} {@var{NewMdl} =} step (@var{mdl})
    ## @deftypefnx {LinearModel} {@var{NewMdl} =} step (@var{mdl}, @var{Name}, @var{Value})
    ##
    ## Improve a fitted linear regression model by one or more steps of
    ## stepwise term selection.
    ##
    ## @code{step} examines whether adding or removing a single term from
    ## @var{mdl} improves the fit, and returns the resulting model as
    ## @var{NewMdl}.  The original model @var{mdl} is never modified.  Unlike
    ## @code{stepwiselm}, @code{step} performs only one such improvement step
    ## by default; pass @code{'NSteps'} to allow more.
    ##
    ## @code{step} accepts the same @code{'Criterion'}, @code{'PEnter'},
    ## @code{'PRemove'}, @code{'NSteps'}, @code{'Verbose'}, @code{'Lower'},
    ## and @code{'Upper'} Name-Value options as @code{stepwiselm}, with the
    ## same defaults, except @code{'NSteps'} defaults to @code{1} rather than
    ## unlimited.  @var{mdl}'s own predictors, weights, excluded observations,
    ## and categorical variable settings are carried over automatically as
    ## the starting point for the search.
    ##
    ## @code{step} is not available for a model fitted with robust
    ## regression.
    ##
    ## @seealso{stepwiselm, addTerms, removeTerms}
    ## @end deftypefn
    function NewMdl = step (mdl, varargin)
      if (nargin < 1)
        error ("step: Not enough input arguments.");
      endif
      if (! isempty (mdl.Robust))
        error ("step: The STEP method is not available with a robust fit.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error ("step: Name-Value arguments must be in pairs.");
      endif

      has_nsteps = false;
      for i = 1:2:numel (varargin) - 1
        if ((ischar (varargin{i}) || isstring (varargin{i})) ...
            && strcmpi (char (varargin{i}), 'NSteps'))
          has_nsteps = true;
        endif
      endfor
      extra = varargin;
      if (! has_nsteps)
        extra = [extra, {'NSteps', 1}];
      endif

      cat_vars = {};
      if (! isempty (mdl.CatLevelInfo) && isfield (mdl.CatLevelInfo, 'names'))
        cat_vars = mdl.CatLevelInfo.names;
      endif

      nv_list = {'PredictorVars', mdl.PredictorNames};
      if (! isempty (mdl.OrigOpts.Weights))
        nv_list = [nv_list, {'Weights', mdl.OrigOpts.Weights}];
      endif
      if (! isempty (mdl.OrigOpts.Exclude))
        nv_list = [nv_list, {'Exclude', mdl.OrigOpts.Exclude}];
      endif
      if (! isempty (cat_vars))
        nv_list = [nv_list, {'CategoricalVars', cat_vars}];
      endif

      if (mdl.HasIntercept)
        formula_str = [mdl.ResponseName, ' ~ ', mdl.Formula.LinearPredictor];
      elseif (isempty (mdl.Formula.LinearPredictor))
        formula_str = [mdl.ResponseName, ' ~ -1'];
      else
        formula_str = [mdl.ResponseName, ' ~ ', mdl.Formula.LinearPredictor, ' - 1'];
      endif

      NewMdl = stepwiselm (mdl.Variables, formula_str, nv_list{:}, extra{:});
    endfunction

  endmethods

  methods (Access = public, Static, Hidden)

    ## weighted least-squares via pivoted QR; returns fit struct
    function fit = lm_fit (X, y, w, compute_H)
      if (nargin < 4)
        compute_H = true;
      endif
      n = rows (X);
      p = columns (X);
      w = w(:);

      W_sqrt = sqrt (w);
      Xw     = X .* W_sqrt;
      yw     = y .* W_sqrt;

      [Q, R, Pperm] = qr (Xw, 0);

      if (isvector (Pperm))
        P_vec = double (Pperm(:)');
      else
        [~, P_vec] = max (Pperm, [], 1);
      endif

      dr = abs (diag (R));
      if (isempty (dr) || dr(1) == 0)
        rank_X = 0;
      else
        tol    = max (size (Xw)) * eps (dr(1));
        rank_X = sum (dr > tol);
      endif

      beta        = zeros (p, 1);
      active_cols = P_vec(1:rank_X);

      if (rank_X > 0)
        R11   = R(1:rank_X, 1:rank_X);
        Q1    = Q(:, 1:rank_X);
        gamma = R11 \ (Q1' * yw);
        beta(active_cols) = gamma;
      endif

      Fitted = X * beta;
      Raw    = y - Fitted;

      n_eff   = sum (w > 0);
      SumLogW = sum (log (w(w > 0)));
      SSE     = sum (w .* Raw.^2);
      wmean   = sum (w .* y) / max (sum (w), eps);
      SST     = sum (w .* (y - wmean).^2);
      SSR     = SST - SSE;

      DFE = n_eff - rank_X;
      if (DFE > 0)
        MSE  = SSE / DFE;
        RMSE = sqrt (MSE);
      else
        MSE  = NaN;
        RMSE = NaN;
      endif

      CovBeta = zeros (p, p);
      if (rank_X > 0)
        R11_inv = R11 \ eye (rank_X);
        CovBeta(active_cols, active_cols) = MSE * (R11_inv * R11_inv');
      endif

      ## Compute hat matrix and leverage
      if (rank_X > 0)
        leverage = sum (Q1.^2, 2);
        if (compute_H)
          Q1t = Q1 * Q1';
          H   = (Q1t ./ W_sqrt) .* W_sqrt';
        else
          H = [];
        endif
      else
        leverage = zeros (n, 1);
        if (compute_H)
          H = zeros (n, n);
        else
          H = [];
        endif
      endif

      fit.beta        = beta;
      fit.H           = H;
      fit.leverage    = leverage;
      fit.SSE         = SSE;
      fit.SSR         = SSR;
      fit.SST         = SST;
      fit.DFE         = DFE;
      fit.MSE         = MSE;
      fit.RMSE        = RMSE;
      fit.CovBeta     = CovBeta;
      fit.rank_X      = rank_X;
      fit.active_cols = active_cols;
      fit.Fitted      = Fitted;
      fit.Raw         = Raw;
      fit.SumLogW     = SumLogW;
    endfunction

    function crit = lm_criteria (fit, n_obs, has_intercept)
      p   = fit.rank_X;
      SSE = fit.SSE;
      SSR = fit.SSR;
      SST = fit.SST;
      DFE = fit.DFE;
      MSE = fit.MSE;

      LogLikelihood = -(n_obs / 2) * (1 + log (2 * pi * SSE / n_obs)) + 0.5 * fit.SumLogW;

      AIC  = -2 * LogLikelihood + 2 * p;
      dAIC = n_obs - p - 1;
      if (dAIC > 0)
        AICc = AIC + (2 * p * (p + 1)) / dAIC;
      else
        AICc = Inf;
      endif
      BIC  = -2 * LogLikelihood + p * log (n_obs);
      CAIC = BIC + p;

      R2_ord = SSR / max (SST, eps);
      if (n_obs > 1 && DFE > 0)
        R2_adj = 1 - (SSE / DFE) / (SST / (n_obs - 1));
      else
        R2_adj = NaN;
      endif

      if (has_intercept && p > 1)
        df1   = p - 1;
        Fstat = (SSR / df1) / max (MSE, eps);
      elseif (! has_intercept && p > 0)
        df1   = p;
        Fstat = (SSR / df1) / max (MSE, eps);
      else
        df1   = 0;
        Fstat = NaN;
      endif

      if (df1 > 0 && DFE > 0 && Fstat >= 0)
        Fpval = betainc (DFE / (DFE + df1 * Fstat), DFE / 2, df1 / 2);
      else
        Fpval = NaN;
      endif

      crit.LogLikelihood = LogLikelihood;
      crit.AIC           = AIC;
      crit.AICc          = AICc;
      crit.BIC           = BIC;
      crit.CAIC          = CAIC;
      crit.Rsquared      = R2_ord;
      crit.AdjRsquared   = R2_adj;
      crit.Fstat         = Fstat;
      crit.Fpval         = Fpval;
    endfunction

    function info = sw_extract (mdl0)
      pred_names = mdl0.PredictorNames;
      p_raw      = mdl0.NumPredictors;
      cat_info   = mdl0.CatLevelInfo;
      enc_names  = mdl0.EncPredictorNames;

      y_sub = mdl0.ResponseVector (mdl0.SubsetMask);
      w_sub = mdl0.WeightVector (mdl0.SubsetMask);

      if (! isempty (mdl0.EncodedPredMatrix))
        X_enc_sub = mdl0.EncodedPredMatrix;
      else
        tbl_sub = mdl0.Variables (mdl0.SubsetMask, :);
        X_raw   = zeros (rows (tbl_sub), p_raw);
        for j = 1:p_raw
          col = tbl_sub.(pred_names{j});
          if (iscell (col))
            ci       = find (strcmp (cat_info.names, pred_names{j}));
            levels_j = cat_info.levels{ci};
            codes    = zeros (rows (tbl_sub), 1);
            for k = 1:numel (levels_j)
              codes(strcmp (col, levels_j{k})) = k;
            endfor
            X_raw(:, j) = codes;
          elseif (isa (col, 'categorical'))
            ci       = find (strcmp (cat_info.names, pred_names{j}));
            levels_j = cat_info.levels{ci};
            [~, X_raw(:, j)] = ismember (cellstr (col), levels_j);
          else
            X_raw(:, j) = double (col);
          endif
        endfor
        X_enc_sub = reencode_predictors (X_raw, pred_names, cat_info, enc_names);
      endif

      info.terms_enc     = mdl0.TermsMatrix;
      info.cat_info      = cat_info;
      info.enc_names     = enc_names;
      info.pred_names    = pred_names;
      info.p_raw         = p_raw;
      info.X_enc         = X_enc_sub;
      info.y             = y_sub;
      info.w             = w_sub;
      info.has_intercept = mdl0.HasIntercept;
      info.n_obs         = mdl0.NumObservations;
      info.orig_opts     = mdl0.OrigOpts;
      info.variables     = mdl0.Variables;
      info.response_name = mdl0.ResponseName;
    endfunction

  endmethods

endclassdef

function opts = lm_parse_nv (nv_args)

  opt_names = {'Intercept', 'Weights', 'Exclude', 'RobustOpts', ...
               'VarNames', 'CategoricalVars', 'ResponseVar', 'PredictorVars'};
  def_vals  = {true, [], [], [], {}, [], '', {}};
  [intercept, weights, exclude, robustopts, varnames, catvars, ...
   respvar, predvars, rem_args] = parsePairedArguments (opt_names, def_vals, nv_args);

  opts.PredictorVarsGiven = false;
  for i = 1:2:numel (nv_args)
    if ((ischar (nv_args{i}) || isstring (nv_args{i})) ...
        && strcmpi (char (nv_args{i}), 'PredictorVars'))
      opts.PredictorVarsGiven = true;
    endif
  endfor

  if (! isempty (rem_args))
    error ("LinearModel: Unknown option '%s'.", rem_args{1});
  endif

  opts.Intercept       = logical (intercept);
  opts.Exclude         = exclude;
  opts.CategoricalVars = catvars;
  opts.ResponseVar     = char (respvar);

  if (isempty (robustopts) || (ischar (robustopts) && strcmpi (robustopts, 'off')) ...
      || (islogical (robustopts) && ! robustopts))
    opts.RobustOpts = [];
  else
    rname = 'bisquare';
    rtune = [];
    if (isstruct (robustopts))
      if (isfield (robustopts, 'RobustWgtFun') && ! isempty (robustopts.RobustWgtFun))
        rname = robustopts.RobustWgtFun;
      endif
      if (isfield (robustopts, 'Tune'))
        rtune = robustopts.Tune;
      endif
    elseif (is_function_handle (robustopts))
      rname = robustopts;
    elseif (ischar (robustopts) && ! strcmpi (robustopts, 'on'))
      rname = robustopts;
    elseif (! (ischar (robustopts) && strcmpi (robustopts, 'on')))
      error ("LinearModel: invalid RobustOpts value.");
    endif

    if (is_function_handle (rname))
      wfun = rname;
      if (isempty (rtune))
        rtune = 1;
      endif
    else
      switch (lower (char (rname)))
        case 'andrews'
          wfun = @(r) (abs (r) < pi) .* sin (max (sqrt (eps), abs (r))) ...
                      ./ max (sqrt (eps), abs (r));
          def_tune = 1.339;
        case 'bisquare'
          wfun = @(r) (abs (r) < 1) .* (1 - r.^2).^2;
          def_tune = 4.685;
        case 'cauchy'
          wfun = @(r) 1 ./ (1 + r.^2);
          def_tune = 2.385;
        case 'fair'
          wfun = @(r) 1 ./ (1 + abs (r));
          def_tune = 1.400;
        case 'huber'
          wfun = @(r) 1 ./ max (1, abs (r));
          def_tune = 1.345;
        case 'logistic'
          wfun = @(r) tanh (max (sqrt (eps), abs (r))) ./ max (sqrt (eps), abs (r));
          def_tune = 1.205;
        case 'ols'
          wfun = @(r) ones (size (r));
          def_tune = 1;
        case 'talwar'
          wfun = @(r) 1 * (abs (r) < 1);
          def_tune = 2.795;
        case 'welsch'
          wfun = @(r) exp (-(r.^2));
          def_tune = 2.985;
        otherwise
          error ("LinearModel: unrecognised RobustWgtFun '%s'.", char (rname));
      endswitch
      if (isempty (rtune))
        rtune = def_tune;
      endif
    endif
    opts.RobustOpts = struct ('WgtFun', wfun, 'Tune', rtune);
  endif

  if (isempty (weights))
    opts.Weights = [];
  else
    opts.Weights = double (weights(:));
  endif

  if (isempty (varnames))
    opts.VarNames = {};
  else
    opts.VarNames = cellstr (varnames);
  endif

  if (isempty (predvars))
    opts.PredictorVars = {};
  elseif (ischar (predvars) || iscellstr (predvars) || isstring (predvars))
    opts.PredictorVars = cellstr (predvars);
  else
    opts.PredictorVars = predvars;
  endif

endfunction

## observation-level influence statistics; returns D struct
function D = lm_diagnostics (X, y, fit, w)
  n    = rows (X);
  p    = fit.rank_X;
  h    = fit.leverage;
  Raw  = fit.Raw;
  DFE  = fit.DFE;
  MSE  = fit.MSE;
  RMSE = fit.RMSE;

  S2_i = (DFE * MSE - w .* Raw.^2 ./ max (1 - h, eps)) / max (DFE - 1, 1);

  r_std = Raw ./ max (RMSE .* sqrt (max (1 - h, eps)), eps);
  r_stu = Raw ./ max (sqrt (max (S2_i, eps)) .* sqrt (max (1 - h, eps)), eps);

  CooksDistance = (w / max (p, 1)) .* r_std.^2 .* h ./ max (1 - h, eps);
  Dffits        = r_stu .* sqrt (h ./ max (1 - h, eps)) .* sqrt (w);
  CovRatio      = (S2_i ./ max (MSE, eps)).^p ./ max (1 - h, eps);

  p_full   = columns (X);  
  active   = fit.active_cols;
  CovB_act = fit.CovBeta(active, active);
  XtXinv_d = diag (CovB_act) / max (MSE, eps);
  Dfbetas  = NaN (n, p_full);

  if (p > 0)
    for i = 1:n
      xi_act     = X(i, active)';
      infl       = (CovB_act / max (MSE, eps)) * xi_act;
      denom_base = (1 - h(i)) * sqrt (max (S2_i(i), eps));
      for jj = 1:p
        se_jj = sqrt (max (XtXinv_d(jj), eps));
        Dfbetas(i, active(jj)) = infl(jj) * Raw(i) / max (denom_base * se_jj, eps);
      endfor
    endfor
  endif

  D.Leverage      = h;
  D.CooksDistance = CooksDistance;
  D.Dffits        = Dffits;
  D.S2_i          = S2_i;
  D.CovRatio      = CovRatio;
  D.Dfbetas       = Dfbetas;
  D.HatMatrix     = fit.H;
endfunction

function mdl2 = lm_refit (mdl, new_terms)
  opts    = mdl.OrigOpts;
  has_int = any (all (new_terms(:, 1:end-1) == 0, 2));

  if (! isempty (mdl.CatLevelInfo) && isfield (mdl.CatLevelInfo, 'names') ...
      && ! isempty (mdl.CatLevelInfo.names))
    cat_vars = mdl.CatLevelInfo.names;
  else
    cat_vars = opts.CategoricalVars;
  endif

  nv_list = {'Intercept', has_int};
  if (! isempty (opts.Weights))
    nv_list = [nv_list, {'Weights', opts.Weights}];
  endif
  if (! isempty (opts.Exclude))
    nv_list = [nv_list, {'Exclude', opts.Exclude}];
  endif
  if (! isempty (cat_vars))
    nv_list = [nv_list, {'CategoricalVars', cat_vars}];
  endif

  mdl2 = fitlm (mdl.Variables, mdl.ResponseName, new_terms, nv_list{:});
endfunction

function [ax, mdl, args] = lm_plot_axes (this, rest)

  if (isscalar (this) && isgraphics (this, 'axes'))
    ax   = this;
    mdl  = rest{1};
    args = rest(2:end);
  else
    ax   = [];
    mdl  = this;
    args = rest;
  endif
endfunction

function props = lm_plot_props (nv_args)

  opt_names = {'Color', 'Marker', 'MarkerSize', 'MarkerEdgeColor', ...
               'MarkerFaceColor', 'LineWidth'};
  def_vals  = {[0.1490, 0.5490, 0.8660], 'x', 6, 'auto', 'none', 0.5};
  [color, marker, markersize, mec, mfc, lw, rem_args] = ...
    parsePairedArguments (opt_names, def_vals, nv_args);
  if (! isempty (rem_args))
    error ("lm_plot_props: unrecognized property '%s'.", rem_args{1});
  endif
  props.Color           = color;
  props.Marker          = marker;
  props.MarkerSize      = markersize;
  props.MarkerEdgeColor = mec;
  props.MarkerFaceColor = mfc;
  props.LineWidth       = lw;
endfunction

function h = lm_plot_data (ax, xdata, ydata, props)
  h = plot (ax, xdata, ydata, ...
            'LineStyle',      'none', ...
            'Color',          props.Color, ...
            'Marker',         props.Marker, ...
            'MarkerSize',     props.MarkerSize, ...
            'MarkerEdgeColor', props.MarkerEdgeColor, ...
            'MarkerFaceColor', props.MarkerFaceColor, ...
            'LineWidth',      props.LineWidth);
endfunction

function loc = lm_legend_corner (xdata, ydata)
  xr   = xdata(! isnan (xdata));
  yr   = ydata(! isnan (ydata));
  xmid = (min (xr) + max (xr)) / 2;
  ymid = (min (yr) + max (yr)) / 2;
  counts = [sum(xr >= xmid & yr >= ymid), sum(xr < xmid & yr >= ymid), ...
            sum(xr >= xmid & yr < ymid),  sum(xr < xmid & yr < ymid)];
  locs = {'northeast', 'northwest', 'southeast', 'southwest'};
  [~, best_idx] = min (counts);
  loc = locs{best_idx};
endfunction

function c_row = lm_interaction_row (X_act, fix_cols, fix_vals, pred, cinfo, ename, terms)
  X_rows = X_act;
  for f = 1:numel (fix_cols)
    X_rows(:, fix_cols(f)) = fix_vals(f);
  endfor
  X_enc = reencode_predictors (X_rows, pred, cinfo, ename);
  D     = build_design (terms, X_enc);
  c_row = mean (D, 1);
endfunction

function col_str = lm_col_to_str (col)
  if (iscell (col))
    col_str = col;
  elseif (isa (col, 'categorical'))
    col_str = cellstr (col);
  elseif (isa (col, 'string'))
    col_str = cellstr (col);
  else
    col_str = cellstr (num2str (col(:)));
  endif
endfunction

function [X_act, is_cat, cat_lvls] = lm_encode_active_predictors (mdl, act, pred, cinfo)
  n_act    = sum (act);
  p        = numel (pred);
  X_act    = zeros (n_act, p);
  is_cat   = false (1, p);
  cat_lvls = cell (1, p);

  for k = 1:p
    ci = [];
    if (! isempty (cinfo) && isfield (cinfo, 'names') ...
        && ! isempty (cinfo.names))
      ci = find (strcmp (cinfo.names, pred{k}));
    endif
    col = mdl.Variables{act, pred{k}};
    if (! isempty (ci))
      is_cat(k)   = true;
      levels_k    = cinfo.levels{ci};
      cat_lvls{k} = levels_k;
      col_str = lm_col_to_str (col);
      codes = zeros (n_act, 1);
      for L = 1:numel (levels_k)
        codes(strcmp (col_str, char (levels_k{L}))) = L;
      endfor
      X_act(:,k) = codes;
    else
      X_act(:,k) = double (col(:));
    endif
  endfor
endfunction

function used = lm_predictors_in_model (pred_all, cinfo, ename)
  n    = numel (pred_all);
  used = false (1, n);
  for k = 1:n
    ci = [];
    if (! isempty (cinfo) && ! isempty (cinfo.names))
      ci = find (strcmp (cinfo.names, pred_all{k}));
    endif
    if (isempty (ci))
      used(k) = any (strcmp (ename, pred_all{k}));
    else
      levels_k = cinfo.levels{ci};
      for L = 2:numel (levels_k)
        nm = sprintf ("%s_%s", pred_all{k}, char (levels_k{L}));
        if (any (strcmp (ename, nm)))
          used(k) = true;
        endif
      endfor
    endif
  endfor
endfunction

function C = lm_effects_contrasts (mdl)
  if (! any (any (mdl.TermsMatrix(:, 1:end-1) != 0)))
    C = [];
    return;
  endif

  pred_all = mdl.PredictorNames;
  cinfo    = mdl.CatLevelInfo;
  ename    = mdl.EncPredictorNames;
  pred     = pred_all (lm_predictors_in_model (pred_all, cinfo, ename));
  p        = numel (pred);
  act      = mdl.ObservationInfo.Subset;

  [X_act, is_cat, cat_lvls] = lm_encode_active_predictors (mdl, act, pred, cinfo);

  C = zeros (p, mdl.NumCoefficients);

  for j = 1:p
    if (is_cat(j))
      x_lo = 1;
      x_hi = numel (cat_lvls{j});
    else
      x_lo = min (X_act(:,j));
      x_hi = max (X_act(:,j));
    endif

    X_hi_rows = X_act;  X_hi_rows(:,j) = x_hi;
    X_lo_rows = X_act;  X_lo_rows(:,j) = x_lo;

    X_hi_enc = reencode_predictors (X_hi_rows, pred, cinfo, ename);
    X_lo_enc = reencode_predictors (X_lo_rows, pred, cinfo, ename);
    D_hi     = build_design (mdl.TermsMatrix, X_hi_enc);
    D_lo     = build_design (mdl.TermsMatrix, X_lo_enc);

    C(j,:) = mean (D_hi - D_lo, 1);
  endfor
endfunction

function IC = lm_interaction_contrasts (mdl)
  pred_all = mdl.PredictorNames;
  cinfo    = mdl.CatLevelInfo;
  ename    = mdl.EncPredictorNames;
  pred     = pred_all (lm_predictors_in_model (pred_all, cinfo, ename));
  p        = numel (pred);

  IC = struct ('OwnGridRows', {cell(1, p)}, 'Pairs', {cell(p, p)});
  if (p < 2)
    return;
  endif

  act   = mdl.ObservationInfo.Subset;
  terms = mdl.TermsMatrix;
  tpred = terms(:, 1:p);

  [X_act, is_cat, cat_lvls] = lm_encode_active_predictors (mdl, act, pred, cinfo);

  grids = cell (1, p);
  for j = 1:p
    if (is_cat(j))
      n_lv = numel (cat_lvls{j});
      grids{j} = (1:n_lv)';
      rows_j = zeros (n_lv, mdl.NumCoefficients);
      for L = 1:n_lv
        rows_j(L,:) = lm_interaction_row (X_act, j, L, pred, cinfo, ename, terms);
      endfor
      IC.OwnGridRows{j} = rows_j;
    else
      lo = min (X_act(:,j));
      hi = max (X_act(:,j));
      grids{j} = [lo; (lo + hi) / 2; hi];
    endif
  endfor

  for j1 = 1:p
    for j2 = 1:p
      if (j1 == j2)
        continue;
      endif

      shared = tpred(:,j1) > 0 & tpred(:,j2) > 0;
      if (! any (shared))
        continue;
      endif

      grid1 = grids{j1};

      if (is_cat(j2))
        grid2 = grids{j2};
      else
        lo2  = min (X_act(:,j2));
        hi2  = max (X_act(:,j2));
        mid2 = (lo2 + hi2) / 2;
        deg2 = max (tpred(shared, j2));
        extra_n = max (0, deg2 + 1 - 3);
        if (extra_n > 0)
          q = linspace (lo2, hi2, extra_n + 2);
          grid2 = sort ([lo2; mid2; hi2; q(2:end-1)']);
        else
          grid2 = [lo2; mid2; hi2];
        endif
      endif

      n1   = numel (grid1);
      n2   = numel (grid2);
      rows = zeros (n1 * n2, mdl.NumCoefficients);
      idx  = 0;
      for a = 1:n1
        for b = 1:n2
          idx = idx + 1;
          rows(idx,:) = lm_interaction_row (X_act, [j1, j2], [grid1(a), grid2(b)], ...
                                             pred, cinfo, ename, terms);
        endfor
      endfor

      IC.Pairs{j1,j2} = struct ('grid1', grid1, 'grid2', grid2, 'rows', rows);
    endfor
  endfor
endfunction

function fit = lm_robust_fit (X, y, w, wgtfun, tune)

  n       = rows (X);
  p       = columns (X);
  w       = w(:);
  sw      = sqrt (w);
  SumLogW = sum (log (w(w > 0)));

  Xw   = X .* sw;
  yw   = y .* sw;
  beta = Xw \ yw;

  [~, R] = qr (Xw, 0);
  E = Xw / R;
  h = min (0.9999, sum (E.^2, 2));
  adjfactor = 1 ./ sqrt (max (1 - h, eps));

  DFE = n - p;
  if (DFE <= 0)
    fit.beta          = beta;
    fit.H             = zeros (n, n);
    fit.leverage      = zeros (n, 1);
    fit.SSE           = NaN;
    fit.SSR           = NaN;
    fit.SST           = NaN;
    fit.DFE           = DFE;
    fit.MSE           = NaN;
    fit.RMSE          = NaN;
    fit.CovBeta       = NaN (p, p);
    fit.rank_X        = p;
    fit.active_cols   = 1:p;
    fit.Fitted        = X * beta;
    fit.Raw           = y - fit.Fitted;
    fit.RobustWeights = ones (n, 1);
    fit.SumLogW       = SumLogW;
    return;
  endif
  ols_s = norm (y - X * beta) / sqrt (DFE);

  tiny_s  = 1e-6 * std (y);
  tolD    = sqrt (eps);
  iterlim = 50;
  iter    = 0;
  beta0   = zeros (size (beta));
  wts     = ones (n, 1);

  while (iter == 0 || any (abs (beta - beta0) > tolD * max (abs (beta), abs (beta0))))
    iter = iter + 1;
    if (iter > iterlim)
      break;
    endif
    r    = (y - X * beta) ./ sw;
    radj = r .* adjfactor;

    rs = sort (abs (radj));
    s  = median (rs(max (1, p):end)) / 0.6745;

    wts   = wgtfun(radj / (max (s, tiny_s) * tune));
    beta0 = beta;

    ww   = sqrt (w .* wts);
    beta = (X .* ww) \ (y .* ww);
  endwhile

  r    = (y - X * beta) ./ sw;
  radj = r .* adjfactor;
  rs    = sort (abs (radj));
  mad_s = median (rs(max (1, p):end)) / 0.6745;

  if (all (wts < tolD | wts > 1 - tolD))
    included = wts > 1 - tolD;
    robust_s = norm (r(included)) / sqrt (max (sum (included) - p, eps));
  else
    st  = max (mad_s, tiny_s) * tune;
    u   = radj / st;
    phi = u .* wgtfun(u);

    delta = 0.0001;
    u1   = u - delta;
    phi0 = u1 .* wgtfun(u1);
    u1   = u + delta;
    phi1 = u1 .* wgtfun(u1);
    dphi = (phi1 - phi0) / (2 * delta);

    m1 = mean (dphi);
    m2 = sum ((1 - h) .* phi.^2) / (n - p);
    K  = 1 + (p / n) * (1 - m1) / m1;
    robust_s = K * sqrt (m2) * st / m1;
  endif

  sigma = max (robust_s, sqrt ((ols_s^2 * p^2 + robust_s^2 * n) / (p^2 + n)));

  RI      = R \ eye (p);
  CovBeta = (RI * RI') * sigma^2;

  ww  = sqrt (w .* wts);
  Xwf = X .* ww;
  [Qf, ~] = qr (Xwf, 0);
  Q1t = Qf * Qf';
  H   = (Q1t ./ ww) .* ww';

  Fitted = X * beta;
  Raw    = y - Fitted;
  ybar   = mean (y);
  SSR    = sum ((Fitted - ybar).^2);
  SSE    = sigma^2 * DFE;
  SST    = SSE + SSR;

  fit.beta          = beta;
  fit.H             = H;
  fit.leverage      = h;
  fit.SSE           = SSE;
  fit.SSR           = SSR;
  fit.SST           = SST;
  fit.DFE           = DFE;
  fit.MSE           = sigma^2;
  fit.RMSE          = sigma;
  fit.CovBeta       = CovBeta;
  fit.rank_X        = p;
  fit.active_cols   = 1:p;
  fit.Fitted        = Fitted;
  fit.Raw           = Raw;
  fit.RobustWeights = wts;
  fit.SumLogW       = SumLogW;

endfunction

function [delta, fit_f, fit_r] = anova_delta_sse (X, y, w, cols_reduced, cols_full, fit_f)
  if (nargin < 6 || isempty (fit_f))
    fit_f = LinearModel.lm_fit (X(:, cols_full), y, w, false);
  endif
  fit_r = LinearModel.lm_fit (X(:, cols_reduced), y, w, false);
  delta = fit_r.SSE - fit_f.SSE;
endfunction

function contain_mx = anova_containment (groups, extended)
  nterm = numel (groups);
  factor_list = cell (nterm, 1);
  for k = 1:nterm
    parts = strsplit (groups(k).Name, ':');
    fl = struct ('var', {}, 'exp', {});
    for f = 1:numel (parts)
      p = strsplit (parts{f}, '^');
      if (numel (p) == 2)
        fl(end+1) = struct ('var', p{1}, 'exp', str2double (p{2}));
      else
        fl(end+1) = struct ('var', p{1}, 'exp', 1);
      endif
    endfor
    factor_list{k} = fl;
  endfor

  contain_mx = false (nterm, nterm);
  for i = 1:nterm
    for j = 1:nterm
      if (i == j)
        continue;
      endif
      fi = factor_list{i};
      fj = factor_list{j};
      ok = true;
      for f = 1:numel (fi)
        match = false;
        for g = 1:numel (fj)
          if (strcmp (fi(f).var, fj(g).var))
            if (extended)
              match = (fj(g).exp >= fi(f).exp);
            else
              match = (fj(g).exp == fi(f).exp);
            endif
            if (match)
              break;
            endif
          endif
        endfor
        if (! match)
          ok = false;
          break;
        endif
      endfor
      contain_mx(i, j) = ok;
    endfor
  endfor
endfunction

function X_raw = anova_decode_raw (tbl_sub, pred_names, cat_info)
  p_raw = numel (pred_names);
  X_raw = zeros (rows (tbl_sub), p_raw);
  for j = 1:p_raw
    col = tbl_sub.(pred_names{j});
    if (iscell (col))
      ci       = find (strcmp (cat_info.names, pred_names{j}));
      levels_j = cat_info.levels{ci};
      codes    = zeros (rows (tbl_sub), 1);
      for k = 1:numel (levels_j)
        codes(strcmp (col, levels_j{k})) = k;
      endfor
      X_raw(:, j) = codes;
    elseif (isa (col, 'categorical'))
      ci       = find (strcmp (cat_info.names, pred_names{j}));
      levels_j = cat_info.levels{ci};
      [~, X_raw(:, j)] = ismember (cellstr (col), levels_j);
    else
      X_raw(:, j) = double (col);
    endif
  endfor
endfunction

function p = f_pvalue (dfe, df, Fstat)
  p = betainc (dfe ./ (dfe + df .* Fstat), dfe / 2, df / 2);
endfunction

%!demo
%!
%! ## Simple linear regression with a single predictor.
%! ## Ten runners record their weekly training distance and their finish
%! ## time in a 10k race. We fit a straight line through this data and
%! ## look at the fitted coefficients, then use predict to estimate the
%! ## finish time for a runner who trains a distance not in the sample.
%! Distance = [10; 15; 20; 25; 30; 35; 40; 45; 50; 55];
%! Time     = [58; 55; 52; 50; 47; 45; 43; 41; 40; 38];
%! X = Distance;
%! y = Time;
%!
%! ## Fit the model and inspect the estimated slope and intercept.
%! mdl = fitlm (X, y)
%!
%! ## Predict the finish time for a runner training 32 km per week.
%! ypred = predict (mdl, 32)

%!demo
%!
%! ## Multiple linear regression with two predictors, followed by a
%! ## confidence interval on the coefficients.
%! ## Thirteen coffee shops report their weekly foot traffic and the
%! ## number of items on their menu, along with weekly revenue. We fit a
%! ## model with both predictors, then use coefCI to see how precisely
%! ## each coefficient is estimated.
%! Traffic = [120; 150; 90; 200; 175; 60; 220; 140; 100; 190; 80; 210; 130];
%! MenuSize = [8; 12; 6; 15; 10; 5; 18; 9; 7; 14; 6; 16; 11];
%! Revenue = [1450; 1820; 1010; 2400; 2050; 700; 2650; 1700; 1150; ...
%!            2300; 900; 2500; 1600];
%! X = [Traffic, MenuSize];
%! y = Revenue;
%!
%! ## Fit the model with both predictors together.
%! mdl = fitlm (X, y)
%!
%! ## Check how tight the 95% confidence interval is on each coefficient.
%! ci = coefCI (mdl)

%!demo
%!
%! ## Growing a model with addTerms and predicting with the richer model.
%! ## We model fuel economy from the carsmall data set using weight and
%! ## horsepower as main effects only. addTerms then brings in the
%! ## weight-horsepower interaction without needing to refit by hand, and
%! ## predict shows how the estimate for a new car changes once that
%! ## interaction is included.
%! load carsmall
%! X = [Weight, Horsepower];
%! y = MPG;
%!
%! ## Fit the additive model first.
%! mdl = fitlm (X, y);
%!
%! ## Add the interaction between weight and horsepower.
%! mdl2 = addTerms (mdl, 'x1:x2');
%!
%! ## Compare predictions from both models for the same new car.
%! Xnew = [3200, 120];
%! ypred1 = predict (mdl, Xnew)
%! ypred2 = predict (mdl2, Xnew)

%!demo
%!
%! ## Simplifying a model with removeTerms and comparing fit quality.
%! ## We fit the full Hald cement model with all four ingredients, then
%! ## use removeTerms to drop the weakest predictor and refit
%! ## automatically. Comparing SSE before and after shows how little
%! ## explanatory power that ingredient was actually contributing.
%! load hald
%! X = ingredients;
%! y = heat;
%!
%! ## Fit the model with all four ingredients.
%! mdl = fitlm (X, y);
%!
%! ## Drop the third ingredient and refit on the same data.
%! mdl2 = removeTerms (mdl, 'x3');
%!
%! ## Compare how much the error sum of squares changed.
%! sse_full    = mdl.SSE
%! sse_reduced = mdl2.SSE

%!demo
%!
%! ## Testing a linear hypothesis and checking residual autocorrelation.
%! ## Twelve patients are given a drug at different doses over different
%! ## treatment durations, and a recovery score is recorded. coefTest
%! ## checks whether the dose and duration coefficients are actually
%! ## equal, and dwtest separately checks whether the residuals still
%! ## carry a leftover pattern the model failed to capture.
%! Dose     = [10; 15; 20; 25; 30; 35; 12; 18; 22; 28; 32; 38];
%! Duration = [5; 7; 9; 11; 13; 15; 6; 8; 10; 12; 14; 16];
%! Recovery = [42; 48; 55; 60; 68; 74; 45; 52; 58; 65; 71; 78];
%! X = [Dose, Duration];
%! y = Recovery;
%! mdl = fitlm (X, y);
%!
%! ## Test H0: the Dose and Duration coefficients are equal.
%! H = [0 1 -1];
%! [p, F, r] = coefTest (mdl, H)
%!
%! ## Check for autocorrelation left over in the residuals.
%! [pdw, dw] = dwtest (mdl)

%!demo
%!
%! ## Checking residuals against fitted values.
%! ## We fit a mileage model on the carsmall data set using weight and
%! ## horsepower, then plot the raw residuals against the fitted values.
%! ## A pattern in this plot, rather than a random scatter, would suggest
%! ## the linear model is missing some curvature in the relationship.
%! load carsmall
%! X = [Weight, Horsepower];
%! y = MPG;
%! mdl = fitlm (X, y);
%! plotResiduals (mdl, 'fitted')

%!demo
%!
%! ## Spotting influential observations with Cook's distance.
%! ## Sixteen houses are matched by size and age to a sale price, but one
%! ## house was sold far above what its size and age would predict. After
%! ## fitting the model, plotDiagnostics with the cookd option highlights
%! ## that single observation as having outsized influence on the fit.
%! Size  = [80; 95; 110; 120; 65; 140; 100; 130; 90; 150; 75; 105; ...
%!          115; 85; 135; 125];
%! Age   = [5; 10; 3; 8; 20; 2; 15; 6; 12; 1; 18; 9; 4; 14; 7; 11];
%! Price = [200; 230; 260; 280; 150; 320; 240; 300; 210; 340; 170; ...
%!          250; 270; 190; 500; 290];
%! X = [Size, Age];
%! y = Price;
%! mdl = fitlm (X, y);
%! plotDiagnostics (mdl, 'cookd')

%!demo
%!
%! ## Comparing the size of each predictor's effect, alongside a
%! ## hypothesis test on the model as a whole.
%! ## We fit a mileage model on the carsmall data set using weight and
%! ## horsepower. plotEffects draws each coefficient's estimate with its
%! ## confidence interval side by side, and coefTest checks whether
%! ## weight's effect is significantly different from horsepower's.
%! load carsmall
%! X = [Weight, Horsepower];
%! y = MPG;
%! mdl = fitlm (X, y);
%!
%! ## Visualize the relative size of each predictor's effect.
%! plotEffects (mdl)
%!
%! ## Test whether the two coefficients differ significantly.
%! H = [0 1 -1];
%! [p, F, r] = coefTest (mdl, H)

%!shared mdl, X, y, n
%! n = 20;
%! X = [1:n; (1:n).^2]' / n;
%! y = X * [3; -1] + 0.2 * sin ((1:n)');
%! mdl = fitlm (X, y);

%!test
%! ## scalar fit-quality
%! assert_equal (mdl.NumObservations,          20);
%! assert_equal (mdl.NumCoefficients,           3);
%! assert_equal (mdl.NumVariables,              3);
%! assert_equal (mdl.NumPredictors,             2);
%! assert_equal (mdl.NumEstimatedCoefficients,  3);
%! assert_equal (mdl.DFE,                      17);
%! assert_equal (mdl.SSE,  0.386545331386823,   1e-9);
%! assert_equal (mdl.SSR,  583.523874670959,    1e-6);
%! assert_equal (mdl.SST,  583.910420002346,    1e-6);
%! assert_equal (mdl.MSE,  0.0227379606698351,  1e-10);
%! assert_equal (mdl.RMSE, 0.150791116017606,   1e-10);
%! assert_equal (mdl.Rsquared.Ordinary, 0.999338005765704, 1e-10);
%! assert_equal (mdl.Rsquared.Adjusted, 0.999260124091081, 1e-10);
%! assert_equal (mdl.LogLikelihood, 11.0836133807695, 1e-6);
%! assert_equal (mdl.ModelCriterion.AIC,  -16.1672267615389, 1e-6);
%! assert_equal (mdl.ModelCriterion.AICc, -14.6672267615389, 1e-6);
%! assert_equal (mdl.ModelCriterion.BIC,  -13.180029940877,  1e-6);
%! assert_equal (mdl.ModelCriterion.CAIC, -10.180029940877,  1e-6);
%! assert_equal (mdl.ModelFitVsNullModel.Fstat, 12831.4909842738, 1e-4);
%! assert_equal (strcmp (mdl.ModelFitVsNullModel.NullModel, 'constant'), true);

%!test
%! ## constant-only model: SSR is exactly zero, SSE equals SST
%! mc = fitlm (X, y, 'constant');
%! assert_equal (mc.SSE, 583.910420002346, 1e-6);
%! assert_equal (mc.SSR, 0, 1e-12);
%! assert_equal (mc.SSE, mc.SST, 1e-12);

%!test
%! ## coefficient estimates, SE, tStat, names, covariance, schema
%! assert_equal (mdl.Coefficients.Estimate, [0.1161886778; 2.508451491; -0.9788353298], 1e-7);
%! assert_equal (mdl.Coefficients.SE,       [0.112185831;  0.4920818186; 0.02276108523], 1e-8);
%! assert_equal (mdl.Coefficients.tStat,    [1.035680502;  5.097630913; -43.00477415],   1e-6);
%! assert_equal (all (mdl.Coefficients.pValue >= 0 & mdl.Coefficients.pValue <= 1), true);
%! assert_equal (isequal (mdl.CoefficientNames, {'(Intercept)', 'x1', 'x2'}), true);
%! assert_equal (isequal (mdl.CoefficientNames, mdl.Coefficients.Properties.RowNames(:)'), true);
%! assert_equal (size (mdl.CoefficientCovariance), [3, 3]);
%! assert_equal (diag (mdl.CoefficientCovariance), [0.0125857; 0.242145; 0.000518067], 1e-6);
%! assert_equal (width (mdl.Coefficients), 4);
%! assert_equal (isequal (mdl.Coefficients.Properties.VariableNames, ...
%!                  {'Estimate','SE','tStat','pValue'}), true);

%!test
%! ## fitted values, predict(), residual columns (obs 1-3), schema
%! assert_equal (mdl.Fitted, y - mdl.Residuals.Raw, 1e-10);
%! yp = predict (mdl, X);
%! assert_equal (size (yp), [20, 1]);
%! assert_equal (yp(1), 0.192669485827491, 1e-10);
%! assert_equal (yp(2), 0.171266760882256, 1e-10);
%! assert_equal (mdl.Residuals.Raw(1:3),          [0.075624711134088; 0.110592724482880; -0.023756501342530], 1e-10);
%! assert_equal (mdl.Residuals.Pearson(1:3),      [0.501519672586403; 0.733416711830473; -0.157545762442370], 1e-9);
%! assert_equal (mdl.Residuals.Standardized(1:3), [0.632246516521578; 0.844226394951239; -0.172381368754725], 1e-8);
%! assert_equal (mdl.Residuals.Studentized(1:3),  [0.620710275056923; 0.836747864205268; -0.167380843634378], 1e-6);
%! assert_equal (width (mdl.Residuals), 4);
%! assert_equal (isequal (mdl.Residuals.Properties.VariableNames, ...
%!                  {'Raw','Pearson','Studentized','Standardized'}), true);

%!test
%! ## diagnostics
%! H = mdl.Diagnostics.HatMatrix;
%! assert_equal (size (H), [20, 20]);
%! assert_equal (H, H', 1e-10);
%! assert_equal (H * H, H, 1e-8);
%! assert_equal (H(1,1), 0.370779220779221, 1e-10);
%! assert_equal (H(1,2), 0.298051948051948, 1e-10);
%! assert_equal (mdl.Diagnostics.Leverage(1:3),      [0.370779220779221; 0.245283663704716; 0.164718614718615], 1e-10);
%! assert_equal (mdl.Diagnostics.CooksDistance(1:3), [0.078517048682575; 0.077211407930332; 0.001953301452841], 1e-8);
%! assert_equal (mdl.Diagnostics.S2_i(1:3),          [0.023591009857798; 0.023146223303229; 0.024116854077430], 1e-8);
%! assert_equal (mdl.Diagnostics.CovRatio(1:3),      [1.774933176573401; 1.397661919176034; 1.428481535363283], 1e-6);
%! assert_equal (mdl.Diagnostics.Dffits(1:3),        [0.476480465355394; 0.477020506700835; -0.074329411030064], 1e-6);
%! assert_equal (size (mdl.Diagnostics.Dfbetas), [20, 3]);
%! assert_equal (width (mdl.Diagnostics), 7);
%! assert_equal (isequal (mdl.Diagnostics.Properties.VariableNames, ...
%!                  {'Leverage','CooksDistance','Dffits','S2_i', ...
%!                   'CovRatio','Dfbetas','HatMatrix'}), true);

%!test
%! ## ObservationInfo, VariableInfo, names, Formula, Variables
%! assert_equal (width (mdl.ObservationInfo), 4);
%! assert_equal (height (mdl.ObservationInfo), 20);
%! assert_equal (isequal (mdl.ObservationInfo.Properties.VariableNames, ...
%!                  {'Weights','Excluded','Missing','Subset'}), true);
%! assert_equal (all (mdl.ObservationInfo.Weights == 1), true);
%! assert_equal (all (mdl.ObservationInfo.Subset == ...
%!              (! mdl.ObservationInfo.Missing & ! mdl.ObservationInfo.Excluded)), true);
%! assert_equal (width  (mdl.VariableInfo), 4);
%! assert_equal (height (mdl.VariableInfo), 3);
%! assert_equal (isequal (mdl.VariableInfo.Properties.VariableNames, ...
%!                  {'Class','Range','InModel','IsCategorical'}), true);
%! assert_equal (mdl.VariableInfo.InModel(strcmp (mdl.VariableNames, 'y')), false);
%! assert_equal (all (mdl.VariableInfo.InModel(! strcmp (mdl.VariableNames, 'y'))), true);
%! assert_equal (mdl.ResponseName, 'y');
%! assert_equal (isequal (mdl.PredictorNames, {'x1','x2'}), true);
%! assert_equal (isequal (mdl.VariableNames, {'x1','x2','y'}), true);
%! assert_equal (mdl.Formula.HasIntercept, true);
%! assert_equal (mdl.Formula.LinearPredictor, '1 + x1 + x2');
%! assert_equal (mdl.Formula.NTerms, 3);
%! assert_equal (strcmp (mdl.Variables.Properties.VariableNames{end}, 'y'), true);

%!test
%! ## NaN in predictor drops the row from the fit
%! X2 = X;  X2(2,1) = NaN;
%! m = fitlm (X2, y);
%! assert_equal (m.NumObservations, 19);
%! assert_equal (m.ObservationInfo.Missing(2), true);
%! assert_equal (m.ObservationInfo.Subset(2),  false);
%! assert_equal (isnan (m.Fitted(2)), true);
%! assert_equal (m.SSE, 0.370339572851658, 1e-9);
%! assert_equal (m.SST, 547.616796178045,  1e-6);
%! assert_equal (m.Coefficients.Estimate, ...
%!         [0.0641300185953764; 2.68263140079657; -0.985345792254554], 1e-7);
%! yp = predict (m, X2);
%! assert_equal (isnan (yp(2)), true);
%! assert_equal (! isnan (yp(1)), true);
%! assert_equal (size (m.Diagnostics.HatMatrix), [20, 20]);
%! assert_equal (m.Diagnostics.Leverage(1), 0.488485648300892, 1e-8);
%! assert_equal (m.Diagnostics.CooksDistance(1), 0.38266162627176, 1e-6);

%!test
%! ## NaN in response drops the row but predict still works normally since X has no NaN
%! y3 = y;  y3(5) = NaN;
%! m = fitlm (X, y3);
%! assert_equal (m.NumObservations, 19);
%! assert_equal (m.ObservationInfo.Missing(5), true);
%! assert_equal (isnan (m.Fitted(5)), true);
%! assert_equal (m.SSE, 0.337042910721425, 1e-9);
%! assert_equal (m.SST, 558.654961265991,  1e-6);
%! assert_equal (m.Coefficients.Estimate, ...
%!         [0.145131680993155; 2.4865383021829; -0.979635081226208], 1e-7);
%! yp = predict (m, X);
%! assert_equal (yp(5), -0.45777759499388, 1e-8);
%! assert_equal (yp(1), 0.22047684204099, 1e-8);
%! assert_equal (size (m.Diagnostics.HatMatrix), [20, 20]);
%! assert_equal (m.Diagnostics.Leverage(1), 0.386399650026734, 1e-8);

%!test
%! ## multiple NaN rows drop all affected observations from the fit
%! X4 = X;  X4([2,8,14],2) = NaN;
%! m = fitlm (X4, y);
%! assert_equal (sum (m.ObservationInfo.Missing), 3);
%! assert_equal (m.NumObservations, 17);
%! assert_equal (m.SSE, 0.261285495635633, 1e-9);
%! assert_equal (m.SSR, 527.635694805749,  1e-6);
%! assert_equal (m.SST, 527.896980301385,  1e-6);
%! assert_equal (m.Coefficients.Estimate, ...
%!         [0.0986395043600395; 2.3735792982821; -0.97106191310122], 1e-7);
%! assert_equal (size (m.Diagnostics.HatMatrix), [20, 20]);
%! assert_equal (sum (m.Diagnostics.Leverage), 3, 1e-8);

%!test
%! ## exclude by index and exclude by logical vector give identical results
%! m  = fitlm (X, y, 'Exclude', [3, 7]);
%! excl = false (n, 1);  excl([3, 7]) = true;
%! m2 = fitlm (X, y, 'Exclude', excl);
%! assert_equal (m.NumObservations, 18);
%! assert_equal (sum (m.ObservationInfo.Excluded), 2);
%! assert_equal (isnan (m.Fitted(3)) && isnan (m.Fitted(7)), true);
%! assert_equal (m.Coefficients.Estimate, m2.Coefficients.Estimate, 1e-12);
%! assert_equal (m.Coefficients.Estimate, ...
%!         [0.118938102486219; 2.43606890944554; -0.974833228191174], 1e-7);
%! ype = predict (m);
%! assert_equal (size (ype), [20, 1]);
%! assert_equal (! isnan (ype(3)) && ! isnan (ype(7)), true);
%! [~, yci] = predict (m);
%! assert_equal (yci(1,1), -0.0283122762458446, 1e-10);
%! assert_equal (yci(1,2),  0.412312049343719,  1e-10);
%! assert_equal (size (m.Diagnostics.HatMatrix), [20, 20]);
%! assert_equal (m.Diagnostics.Leverage(1), 0.437780279893411, 1e-8);
%! assert_equal (m.Diagnostics.CooksDistance(1), 0.110112457355807, 1e-7);

%!test
%! ## NaN and exclude together remove both the missing and the excluded row
%! X6 = X;  X6(1,1) = NaN;
%! m = fitlm (X6, y, 'Exclude', [2]);
%! assert_equal (m.NumObservations, 18);
%! assert_equal (m.ObservationInfo.Missing(1),  true);
%! assert_equal (m.ObservationInfo.Excluded(2), true);
%! assert_equal (m.SSE, 0.342515396265007, 1e-9);
%! assert_equal (m.Coefficients.Estimate, ...
%!         [-0.0735450184226009; 3.17679029176988; -1.0045827469016], 1e-7);

%!test
%! ## weighted least squares produces different SSE and stores the weights
%! w = abs (sin ((1:n)')) + 0.1;
%! m = fitlm (X, y, 'Weights', w);
%! assert_equal (m.SSE, 0.363519720897775, 1e-10);
%! assert_equal (m.ObservationInfo.Weights, w, 1e-15);
%! assert_equal (m.SST, 4.419834786423099e+02, 1e-8);
%! [yp, yci] = predict (m, [0.5 0.25; 1.0 1.0]);
%! assert_equal (yp(1),  1.106748776307639, 1e-10);
%! assert_equal (yp(2),  1.593185531572655, 1e-10);
%! assert_equal (yci(1,1), 0.763985050242272, 1e-10);
%! assert_equal (yci(1,2), 1.449512502373006, 1e-10);
%! assert_equal (m.Diagnostics.Leverage(1), 0.421642939812731, 1e-8);
%! assert_equal (m.Diagnostics.Leverage(2), 0.301314342928707, 1e-8);
%! assert_equal (m.Diagnostics.HatMatrix(1,1), 0.421642939812731, 1e-8);
%! assert_equal (m.Diagnostics.CooksDistance(1), 0.0728569335883748, 1e-7);
%! assert_equal (m.Diagnostics.CovRatio(1), 1.96611264276187, 1e-6);

%!test
%! ## uniform weights scale internals but leave point estimates unchanged
%! m = fitlm (X, y, 'Weights', 2 * ones (n, 1));
%! assert_equal (m.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-10);

%!test
%! ## constant linear and default modelspecs behave as expected
%! m  = fitlm (X, y, 'constant');
%! assert_equal (m.NumCoefficients, 1);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! m2 = fitlm (X, y, 'linear');
%! m3 = fitlm (X, y, []);
%! assert_equal (m2.NumCoefficients, 3);
%! assert_equal (m2.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-12);
%! assert_equal (m3.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-12);

%!test
%! ## purequadratic modelspec produces the expected term count
%! m = fitlm (X, y, 'purequadratic');
%! assert_equal (m.NumCoefficients, 5);

%!test
%! ## interactions modelspec term count and coefficients are verified
%! m = fitlm (X, y, 'interactions');
%! assert_equal (m.NumCoefficients, 4);
%! assert_equal (m.SSE, 0.383859187927621, 1e-9);
%! assert_equal (m.Coefficients.Estimate, ...
%!         [0.157640728038039; 2.08542680311791; -0.929682701072813; -0.031208018255475], 1e-7);

%!test
%! ## quadratic modelspec is rank deficient for this design and drops one coefficient
%! m = fitlm (X, y, 'quadratic');
%! assert_equal (m.NumCoefficients, 6);
%! assert_equal (m.SSE, 0.315784637443501, 1e-9);
%! assert_equal (m.Coefficients.Estimate, ...
%!         [0.447436249544699; -2.44859403731902; -0.0121968798776254; ...
%!          -1.36755100280532; 0; 0.0318176901083297], 1e-7);
%! drop = find (m.Coefficients.SE == 0);
%! assert_equal (numel (drop), 1);
%! assert_equal (isnan (m.Coefficients.tStat(drop)), true);

%!test
%! ## full modelspec with two predictors matches interactions exactly
%! m  = fitlm (X, y, 'full');
%! m2 = fitlm (X, y, 'interactions');
%! assert_equal (m.NumCoefficients, 4);
%! assert_equal (m.Coefficients.Estimate, m2.Coefficients.Estimate, 1e-10);
%! assert_equal (m.Coefficients.Estimate, ...
%!         [0.157640728038039; 2.08542680311791; -0.929682701072813; -0.031208018255475], 1e-7);

%!test
%! ## full modelspec with three predictors includes the three way interaction term
%! X3 = [X, cos((1:n)' * pi / n)];
%! m  = fitlm (X3, y, 'full');
%! assert_equal (m.NumCoefficients, 8);
%! assert_equal (any (strcmp (m.CoefficientNames, 'x1:x2:x3')), true);
%! assert_equal (m.SSE, 0.231331066631196, 1e-8);
%! idx3 = find (strcmp (m.CoefficientNames, 'x1:x2:x3'));
%! assert_equal (m.Coefficients.Estimate(idx3), 0.514890561912964, 1e-6);

%!test
%! ## full modelspec without an intercept drops the intercept coefficient
%! m = fitlm (X, y, 'full', 'Intercept', false);
%! assert_equal (m.NumCoefficients, 3);
%! assert_equal (! any (strcmp (m.CoefficientNames, '(Intercept)')), true);
%! assert_equal (m.Coefficients.Estimate, ...
%!         [3.232987312533958; -1.041484635851565; 0.0324190990982863], 1e-7);

%!test
%! ## a p column terms matrix produces a model with no intercept
%! m = fitlm (X, y, [1 0; 0 1]);
%! assert_equal (m.NumCoefficients, 2);
%! assert_equal (! any (strcmp (m.CoefficientNames, '(Intercept)')), true);
%! assert_equal (m.Coefficients.Estimate, [2.96142161317611; -0.997248749443286], 1e-7);

%!test
%! ## a p plus one column terms matrix produces a model with an intercept
%! m = fitlm (X, y, [0 0 0; 1 0 0; 0 1 0]);
%! assert_equal (m.NumCoefficients, 3);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.Coefficients.Estimate, ...
%!         [0.116188677790207; 2.50845149057086; -0.978835329825186], 1e-7);

%!test
%! ## a table with a Wilkinson formula fits the same model and predicts on a table
%! T = table (X(:,1), X(:,2), y, 'VariableNames', {'a','b','resp'});
%! m = fitlm (T, 'resp ~ a + b');
%! assert_equal (m.NumCoefficients, 3);
%! assert_equal (m.ResponseName, 'resp');
%! assert_equal (m.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-8);
%! Xt = table ([0.5;1.0], [0.25;1.0], 'VariableNames', {'a','b'});
%! yp = predict (m, Xt);
%! assert_equal (yp(1), 1.125705590619342, 1e-10);
%! assert_equal (yp(2), 1.645804838535884, 1e-10);

%!test
%! ## a matrix with a Wilkinson formula string fits the same model as the matrix alone
%! m = fitlm (X, y, 'y ~ x1 + x2');
%! assert_equal (m.NumCoefficients, 3);
%! assert_equal (m.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-8);

%!test
%! ## a pure interaction formula keeps EncPredictorNames aligned with TermsMatrix
%! mi = fitlm (X, y, 'y ~ x1:x2');
%! assert_equal (numel (mi.EncPredictorNames), columns (mi.TermsMatrix) - 1);
%! assert_equal (mi.NumCoefficients, 2);
%! assert_equal (mi.CoefficientNames, {'(Intercept)', 'x1:x2'});
%! assert_equal (mi.Coefficients.Estimate, ...
%!   [-0.755813941484483; -0.876953077491396], 1e-9);
%! assert_equal (predict (mi), mi.Fitted, 1e-10);
%! fig = figure ('visible', 'off');
%! h = plot (mi);
%! assert_equal (numel (h), 3);
%! assert_equal (get (get (gca, 'Title'),  'String'), 'Added Variable Plot for x1:x2');
%! assert_equal (get (get (gca, 'XLabel'), 'String'), 'Adjusted x1:x2');
%! assert_equal (get (h(2), 'DisplayName'), 'Fit: y = -0.876953*x');
%! close (fig);

%!test
%! ## a table input with the default formula fits the same model as the matrix
%! T3 = table (X(:,1), X(:,2), y, 'VariableNames', {'x1','x2','y'});
%! m = fitlm (T3);
%! assert_equal (m.ResponseName, 'y');
%! assert_equal (m.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-8);

%!test
%! ## VarNames sets custom names and ResponseVar overrides the response name
%! m  = fitlm (X, y, 'VarNames', {'alpha','beta','resp'});
%! assert_equal (m.ResponseName, 'resp');
%! assert_equal (isequal (m.PredictorNames, {'alpha','beta'}), true);
%! assert_equal (any (strcmp (m.CoefficientNames, 'alpha')), true);
%! assert_equal (any (strcmp (m.CoefficientNames, 'beta')), true);
%! m2 = fitlm (X, y, 'VarNames', {'a','b','r'}, 'ResponseVar', 'r');
%! assert_equal (m2.ResponseName, 'r');

%!test
%! ## a rank deficient design matrix leaves the dropped coefficients as NaN across the board
%! X_rd = [ones(n,1), X, X(:,1)+X(:,2)];
%! m = fitlm (X_rd, y);
%! assert_equal (m.NumCoefficients, 5);
%! assert_equal (m.NumEstimatedCoefficients, 3);
%! drop = find (m.Coefficients.SE == 0);
%! assert_equal (numel (drop), 2);
%! assert_equal (all (isnan (m.Coefficients.tStat(drop))), true);
%! assert_equal (all (isnan (m.Coefficients.pValue(drop))), true);
%! assert_equal (m.SST, 5.839104200023459e+02, 1e-8);
%! assert_equal (all (all (m.CoefficientCovariance(drop,:) == 0)), true);
%! yp = predict (m, X_rd);
%! assert_equal (size (yp), [n, 1]);
%! assert_equal (! any (isnan (yp)), true);
%! assert_equal (size (m.Diagnostics.Dfbetas), [20, 5]);
%! assert_equal (all (isnan (m.Diagnostics.Dfbetas(:, drop)(:))), true);
%! assert_equal (m.Diagnostics.Leverage(1), 0.370779220779221, 1e-8);

%!test
%! ## Intercept=false
%! mni = fitlm (X, y, 'Intercept', false);
%! assert_equal (mni.NumCoefficients, 2);
%! assert_equal (mni.Formula.HasIntercept, false);
%! assert_equal (! any (strcmp (mni.CoefficientNames, '(Intercept)')), true);
%! [yp, yci] = predict (mni, [0.5 0.25; 1.0 1.0]);
%! assert_equal (yp(1), 1.231398619227234, 1e-10);
%! assert_equal (yp(2), 1.964172863732825, 1e-10);
%! assert_equal (yci(1,1), 1.001262470857215, 1e-10);

%!test
%! ## p-column terms matrix
%! m_p = fitlm (X, y, [1 0; 0 1]);
%! assert_equal (m_p.NumCoefficients, 2);
%! assert_equal (! any (strcmp (m_p.CoefficientNames, '(Intercept)')), true);

%!test
%! ## p+1 column terms matrix
%! m_p1 = fitlm (X, y, [0 0 0; 1 0 0; 0 1 0]);
%! assert_equal (m_p1.NumCoefficients, 3);
%! assert_equal (m_p1.CoefficientNames{1}, '(Intercept)');

%!test
%! ## table Wilkinson formula
%! T = table (X(:,1), X(:,2), y, 'VariableNames', {'a','b','resp'});
%! mf = fitlm (T, 'resp ~ a + b');
%! assert_equal (mf.NumCoefficients, 3);
%! assert_equal (mf.ResponseName, 'resp');
%! assert_equal (mf.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-8);
%! Xt = table ([0.5;1.0], [0.25;1.0], 'VariableNames', {'a','b'});
%! yp = predict (mf, Xt);
%! assert_equal (yp(1), 1.125705590619342, 1e-10);
%! assert_equal (yp(2), 1.645804838535884, 1e-10);

%!test
%! ## matrix Wilkinson formula
%! mfm = fitlm (X, y, 'y ~ x1 + x2');
%! assert_equal (mfm.NumCoefficients, 3);
%! assert_equal (mfm.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-8);

%!test
%! ## table default
%! T3 = table (X(:,1), X(:,2), y, 'VariableNames', {'x1','x2','y'});
%! mt = fitlm (T3);
%! assert_equal (mt.ResponseName, 'y');
%! assert_equal (mt.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-8);

%!test
%! ## VarNames sets custom names
%! vn = fitlm (X, y, 'VarNames', {'alpha','beta','resp'});
%! assert_equal (vn.ResponseName, 'resp');
%! assert_equal (isequal (vn.PredictorNames, {'alpha','beta'}), true);
%! assert_equal (any (strcmp (vn.CoefficientNames, 'alpha')), true);
%! assert_equal (any (strcmp (vn.CoefficientNames, 'beta')), true);

%!test
%! ## ResponseVar overrides VarNames
%! rv = fitlm (X, y, 'VarNames', {'a','b','r'}, 'ResponseVar', 'r');
%! assert_equal (rv.ResponseName, 'r');

%!test
%! ## rank-deficient matrix
%! X_rd = [ones(n,1), X, X(:,1)+X(:,2)];
%! m_rd = fitlm (X_rd, y);
%! assert_equal (m_rd.NumCoefficients, 5);
%! assert_equal (m_rd.NumEstimatedCoefficients, 3);
%! drop = find (m_rd.Coefficients.SE == 0);
%! assert_equal (numel (drop), 2);
%! assert_equal (all (isnan (m_rd.Coefficients.tStat(drop))), true);
%! assert_equal (all (isnan (m_rd.Coefficients.pValue(drop))), true);
%! assert_equal (m_rd.SST, 5.839104200023459e+02, 1e-8);
%! assert_equal (all (all (m_rd.CoefficientCovariance(drop,:) == 0)), true);
%! yp = predict (m_rd, X_rd);
%! assert_equal (size (yp), [n, 1]);
%! assert_equal (! any (isnan (yp)), true);
%! assert_equal (yp(1:5), [0.192669485827486; 0.171266760882252; ...
%!                   0.0519805029545;  -0.165189287955771; ...
%!                  -0.480242611848561], 1e-10);

%!test
%! ## predict: ypred and default CI at new points
%! [yp, yci] = predict (mdl, [0.5 0.25; 1.0 1.0]);
%! assert_equal (yp(1),    1.125705590619347, 1e-10);
%! assert_equal (yp(2),    1.645804838535894, 1e-10);
%! assert_equal (yci(1,1), 0.810180780547215, 1e-10);
%! assert_equal (yci(1,2), 1.441230400691478, 1e-10);
%! assert_equal (yci(2,1), 0.858229321851723, 1e-10);
%! assert_equal (yci(2,2), 2.433380355220066, 1e-10);

%!test
%! ## predict: observation interval
%! [~, yci] = predict (mdl, [0.5 0.25; 1.0 1.0], 'Prediction', 'observation');
%! assert_equal (yci(1,1), 0.677632064105988, 1e-10);
%! assert_equal (yci(1,2), 1.573779117132706, 1e-10);

%!test
%! ## predict: alpha 0.01
%! [~, yci] = predict (mdl, [0.5 0.25; 1.0 1.0], 'Alpha', 0.01);
%! assert_equal (yci(1,1), 0.692272619570008, 1e-10);
%! assert_equal (yci(1,2), 1.559138561668685, 1e-10);

%!test
%! ## predict: simultaneous CI
%! [~, yci] = predict (mdl, [0.5 0.25; 1.0 1.0], 'Simultaneous', true);
%! assert_equal (yci(1,1), 0.662572505689338, 1e-10);
%! assert_equal (yci(1,2), 1.588838675549355, 1e-10);

%!test
%! ## predict: no Xnew returns all rows including training
%! [yp, yci] = predict (mdl);
%! assert_equal (size (yp),  [20, 1]);
%! assert_equal (size (yci), [20, 2]);
%! assert_equal (yp(1),    0.192669485827490, 1e-10);
%! assert_equal (yp(2),    0.171266760882255, 1e-10);
%! assert_equal (yci(1,1), -0.001052067982566, 1e-10);
%! assert_equal (yci(1,2),  0.386391039637546, 1e-10);

%!test
%! ## predict: NaN predictor propagates to NaN output and CI
%! [yp, yci] = predict (mdl, [0.5 0.25; NaN 1.0; 1.0 1.0]);
%! assert_equal (yp(1), 1.125705590619347, 1e-10);
%! assert_equal (isnan (yp(2)), true);
%! assert_equal (yp(3), 1.645804838535894, 1e-10);
%! assert_equal (isnan (yci(2,1)), true);
%! assert_equal (isnan (yci(2,2)), true);

%!test
%! ## predict: categorical model predictions at group centres
%! Xc    = [1;1;1;2;2;2;3;3;3];
%! yc    = [2.1;2.3;1.9; 4.1;3.9;4.2; 6.3;5.8;6.1];
%! m_cat = fitlm (Xc, yc, 'linear', 'CategoricalVars', 1);
%! [yp, yci] = predict (m_cat, [1;2;3]);
%! assert_equal (yp(1), 2.099999999999998, 1e-10);
%! assert_equal (yp(2), 4.066666666666667, 1e-10);
%! assert_equal (yp(3), 6.066666666666666, 1e-10);
%! assert_equal (yci(1,1), 1.80971256321669, 1e-10);
%! assert_equal (yci(1,2), 2.3902874367833,  1e-10);
%! assert_equal (yci(2,1), 3.77637922988336, 1e-10);
%! assert_equal (yci(2,2), 4.35695410344997, 1e-10);
%! assert_equal (yci(3,1), 5.77637922988336, 1e-10);
%! assert_equal (yci(3,2), 6.35695410344997, 1e-10);

%!test
%! ## predict: interaction model
%! [yp, yci] = predict (fitlm (X, y, 'interactions'), [0.5 0.25; 1.0 1.0]);
%! assert_equal (yp(1),    0.964032452046850, 1e-10);
%! assert_equal (yp(2),    1.282176811827644, 1e-10);
%! assert_equal (yci(1,1), -0.110763003580605, 1e-10);
%! assert_equal (yci(1,2),  2.038827907674306, 1e-10);

%!test
%! ## predict: weighted model, ypred and CI
%! w  = (1:n)' / sum (1:n);
%! mw = fitlm (X, y, 'Weights', w);
%! [yp, yci] = predict (mw, [0.5 0.25; 1.0 1.0]);
%! assert_equal (yp(1),    1.15833357370544, 1e-10);
%! assert_equal (yp(2),    1.74408669002694, 1e-10);
%! assert_equal (yci(1,1), 0.802165170771357, 1e-10);
%! assert_equal (yci(1,2), 1.51450197663953,  1e-10);
%! assert_equal (yci(2,1), 0.69968979253134,  1e-10);
%! assert_equal (yci(2,2), 2.78848358752254,  1e-10);

%!test
%! ## predict: no-intercept model, ypred and CI
%! mni = fitlm (X, y, 'Intercept', false);
%! [yp, yci] = predict (mni, [0.5 0.25; 1.0 1.0]);
%! assert_equal (yp(1),    1.23139861922723, 1e-10);
%! assert_equal (yp(2),    1.96417286373283, 1e-10);
%! assert_equal (yci(1,1), 1.00126247085704, 1e-10);
%! assert_equal (yci(1,2), 1.46153476759743, 1e-10);
%! assert_equal (yci(2,1), 1.51833851162232, 1e-10);
%! assert_equal (yci(2,2), 2.41000721584333, 1e-10);

%!test
%! ## predict: observation interval combined with simultaneous bound
%! [~, yci] = predict (mdl, [0.5 0.25; 1.0 1.0], ...
%!                      'Prediction', 'observation', 'Simultaneous', true);
%! assert_equal (yci(1,1), 0.46801507632267,  1e-10);
%! assert_equal (yci(1,2), 1.78339610491601,  1e-10);
%! assert_equal (yci(2,1), 0.399032373599106, 1e-10);
%! assert_equal (yci(2,2), 2.89257730347266,  1e-10);

%!test
%! ## output is 2x1 double column vector
%! ysim = random (mdl, [0.5, 0.25; 1.0, 1.0]);
%! assert_equal (size (ysim), [2, 1]);
%! assert_equal (class (ysim), 'double');
%! assert_equal (iscolumn (ysim), true);

%!test
%! ## single row input gives 1x1 output
%! assert_equal (size (random (mdl, [0.5, 0.25])), [1, 1]);

%!test
%! ## predict values are exact and noise added is finite
%! ypred = predict (mdl, [0.5, 0.25; 1.0, 1.0]);
%! ysim  = random (mdl, [0.5, 0.25; 1.0, 1.0]);
%! assert_equal (ypred(1), 1.125705590619342, 1e-10);
%! assert_equal (ypred(2), 1.645804838535884, 1e-10);
%! assert_equal (all (isfinite (ysim - ypred)), true);

%!test
%! ## NaN predictor row gives NaN output, other rows stay finite
%! ysim = random (mdl, [0.5, 0.25; NaN, 1.0; 1.0, 1.0]);
%! assert_equal (size (ysim), [3, 1]);
%! assert_equal (isfinite (ysim(1)), true);
%! assert_equal (isnan (ysim(2)), true);
%! assert_equal (isfinite (ysim(3)), true);

%!test
%! ## two sequential calls produce different output
%! ya = random (mdl, [0.5, 0.25]);
%! yb = random (mdl, [0.5, 0.25]);
%! assert_equal (! isequal (ya, yb), true);

%!test
%! ## random: table input, full training data, weighted and no-intercept
%! ## models all give finite output of the expected size
%! Xt = table ([0.5;1.0], [0.25;1.0], 'VariableNames', {'x1','x2'});
%! mw  = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! mni = fitlm (X, y, 'Intercept', false);
%! assert_equal (size (random (mdl, Xt)), [2, 1]);
%! assert_equal (all (isfinite (random (mdl, Xt))), true);
%! assert_equal (size (random (mdl, X)), [20, 1]);
%! assert_equal (sum (isnan (random (mdl, X))), 0);
%! assert_equal (all (isfinite (random (mw,  [0.5, 0.25; 1.0, 1.0]))), true);
%! assert_equal (all (isfinite (random (mni, [0.5, 0.25; 1.0, 1.0]))), true);

%!test
%! yf = feval (mdl, [0.5 0.25; 1.0 1.0; 0.2 0.04]);
%! assert_equal (size (yf), [3, 1]);
%! assert_equal (class (yf), 'double');
%! assert_equal (yf(1), 1.125705590619342, 1e-10);
%! assert_equal (yf(2), 1.645804838535884, 1e-10);
%! assert_equal (yf(3), 0.578725562711373, 1e-10);
%! assert_equal (yf, predict (mdl, [0.5 0.25; 1.0 1.0; 0.2 0.04]), 1e-10);

%!test
%! yf = feval (mdl, [0.5; 1.0; 0.2], [0.25; 1.0; 0.04]);
%! assert_equal (size (yf), [3, 1]);
%! assert_equal (iscolumn (yf), true);
%! assert_equal (yf, predict (mdl, [0.5 0.25; 1.0 1.0; 0.2 0.04]), 1e-10);

%!test
%! yf = feval (mdl, [0.5, 1.0, 0.2], [0.25, 1.0, 0.04]);
%! assert_equal (size (yf), [1, 3]);
%! assert_equal (isrow (yf), true);
%! assert_equal (yf(:), predict (mdl, [0.5 0.25; 1.0 1.0; 0.2 0.04]), 1e-10);

%!test
%! yf = feval (mdl, 0.5, 0.25);
%! assert_equal (size (yf), [1, 1]);
%! assert_equal (yf, 1.125705590619342, 1e-10);
%! assert_equal (yf, predict (mdl, [0.5 0.25]), 1e-10);

%!test
%! yf = feval (mdl, 0.5, [0.1; 0.2; 0.3]);
%! assert_equal (size (yf), [3, 1]);
%! assert_equal (yf(1), 1.272530890093120, 1e-10);
%! assert_equal (yf(2), 1.174647357110602, 1e-10);
%! assert_equal (yf(3), 1.076763824128083, 1e-10);
%! assert_equal (yf, predict (mdl, [0.5 0.1; 0.5 0.2; 0.5 0.3]), 1e-10);

%!test
%! yf = feval (mdl, [0.1; 0.5; 0.9], 0.25);
%! assert_equal (size (yf), [3, 1]);
%! assert_equal (yf(1), 0.122324994390997, 1e-10);
%! assert_equal (yf(2), 1.125705590619342, 1e-10);
%! assert_equal (yf(3), 2.129086186847688, 1e-10);
%! assert_equal (yf, predict (mdl, [0.1 0.25; 0.5 0.25; 0.9 0.25]), 1e-10);

%!test
%! Weight = [2000;2100;2200;2300;2400;2500;2600;2700;2800;2900;3000; ...
%!           3100;3200;3300;3400;3500;3600;3700;3800;3900];
%! Year   = categorical ([70;70;70;70;70;76;76;76;76;76;76;76;82;82; ...
%!                        82;82;82;82;82;82]);
%! MPG    = [30;29;28;27;26;25;24;23;22;21;20;19;18;17;16;15;14;13;12;11];
%! m  = fitlm (table (MPG, Weight, Year), 'MPG ~ Weight + Year');
%! yf = feval (m, [2500;3000], '76');
%! assert_equal (yf(1), 25.000000000000000, 1e-9);
%! assert_equal (yf(2), 20.000000000000004, 1e-9);
%! yf2 = feval (m, [2500;3000], categorical (70));
%! assert_equal (yf2(1), 24.999999999999996, 1e-9);
%! assert_equal (yf2(2), 20.000000000000000, 1e-9);
%! assert_equal (feval (m, 2800, '82'), 21.999999999999996, 1e-9);
%! assert_equal (isnan (feval (m, 2500, '99')), true);

%!test
%! m = fitlm ((1:n)' / n, 2 * (1:n)' / n + 0.1 * sin ((1:n)'));
%! assert_equal (size (feval (m, 0.5)), [1, 1]);
%! assert_equal (size (feval (m, [0.3; 0.5; 0.9])), [3, 1]);
%! assert_equal (feval (m, 0.5), predict (m, 0.5), 1e-10);
%! assert_equal (feval (m, [0.3; 0.5; 0.9]), predict (m, [0.3; 0.5; 0.9]), 1e-10);

%!test
%! T  = table ([0.5; 1.0; 0.2], [0.25; 1.0; 0.04], 'VariableNames', {'x1', 'x2'});
%! yf = feval (mdl, T);
%! assert_equal (size (yf), [3, 1]);
%! assert_equal (yf, predict (mdl, [0.5 0.25; 1.0 1.0; 0.2 0.04]), 1e-10);

%!test
%! yf = feval (mdl, [0.5 0.25; NaN 1.0; 1.0 1.0]);
%! assert_equal (isfinite (yf(1)), true);
%! assert_equal (isnan (yf(2)), true);
%! assert_equal (isfinite (yf(3)), true);

%!test
%! yf = feval (mdl, [0.5; NaN; 1.0], [0.25; 1.0; 1.0]);
%! assert_equal (isnan (yf(2)), true);
%! yf = feval (mdl, [0.5; 1.0; 1.0], [0.25; NaN; 1.0]);
%! assert_equal (isnan (yf(2)), true);

%!test
%! yf = feval (mdl, X);
%! assert_equal (size (yf), [20, 1]);
%! assert_equal (yf, mdl.Fitted, 1e-10);

%!test
%! m  = fitlm (X, y, 'Intercept', false);
%! yf = feval (m, [0.5 0.25; 1.0 1.0]);
%! assert_equal (yf, predict (m, [0.5 0.25; 1.0 1.0]), 1e-10);
%! assert_equal (feval (m, [0.5; 1.0], [0.25; 1.0]), yf, 1e-10);

%!test
%! m  = fitlm (X, y, 'interactions');
%! yf = feval (m, [0.5 0.25; 1.0 1.0]);
%! assert_equal (yf, predict (m, [0.5 0.25; 1.0 1.0]), 1e-10);
%! assert_equal (feval (m, [0.5; 1.0], [0.25; 1.0]), yf, 1e-10);

%!test
%! m  = fitlm ([1;1;1;2;2;2;3;3;3], [2.1;2.3;1.9;4.1;3.9;4.2;6.3;5.8;6.1], ...
%!             'linear', 'CategoricalVars', 1);
%! yf = feval (m, [1; 2; 3]);
%! assert_equal (yf(1), 2.099999999999998, 1e-10);
%! assert_equal (yf(2), 4.066666666666667, 1e-10);
%! assert_equal (yf(3), 6.066666666666666, 1e-10);

%!test
%! ci = coefCI (mdl);
%! assert_equal (size (ci), [3, 2]);
%! assert_equal (class (ci), 'double');
%! assert_equal (all (ci(:,1) < ci(:,2)), true);
%! assert_equal (ci(1,1), -0.120502736154050,  1e-10);
%! assert_equal (ci(1,2),  0.352880091734465, 1e-10);
%! assert_equal (ci(2,1),  1.470249604061007,  1e-10);
%! assert_equal (ci(2,2),  3.546653377080718,  1e-10);
%! assert_equal (ci(3,1), -1.026857022014626,  1e-10);
%! assert_equal (ci(3,2), -0.930813637635746, 1e-10);

%!test
%! ## midpoints equal estimates
%! ci = coefCI (mdl);
%! t  = tinv (0.975, mdl.DFE);
%! assert_equal ((ci(:,1) + ci(:,2)) / 2, mdl.Coefficients.Estimate, 1e-10);
%! assert_equal (ci(:,2) - ci(:,1), 2 * t * mdl.Coefficients.SE, 1e-10);

%!test
%! assert_equal (coefCI (mdl, 0.05), coefCI (mdl));

%!test
%! ci   = coefCI (mdl);
%! ci01 = coefCI (mdl, 0.01);
%! t01  = tinv (0.995, mdl.DFE);
%! assert_equal (size (ci01), [3, 2]);
%! assert_equal (ci01(1,1), -0.208951721610638, 1e-10);
%! assert_equal (ci01(1,2),  0.441329077191052, 1e-10);
%! assert_equal (ci01(2,1),  1.08228494564489,  1e-10);
%! assert_equal (ci01(2,2),  3.934618035496833, 1e-10);
%! assert_equal (ci01(3,1), -1.044802201703589, 1e-10);
%! assert_equal (ci01(3,2), -0.912868457946783, 1e-10);
%! assert_equal (all ((ci01(:,2) - ci01(:,1)) > (ci(:,2) - ci(:,1))), true);
%! assert_equal (ci01(:,2) - ci01(:,1), 2 * t01 * mdl.Coefficients.SE, 1e-10);

%!test
%! ci0 = coefCI (mdl, 0);
%! assert_equal (all (ci0(:,1) == -Inf), true);
%! assert_equal (all (ci0(:,2) == +Inf), true);

%!test
%! ## alpha=1 collapses to point estimates
%! ci1 = coefCI (mdl, 1);
%! assert_equal (ci1(:,1), mdl.Coefficients.Estimate, 1e-10);
%! assert_equal (ci1(:,2), mdl.Coefficients.Estimate, 1e-10);

%!test
%! m  = fitlm (X, y, 'Intercept', false);
%! ci = coefCI (m);
%! t  = tinv (0.975, m.DFE);
%! assert_equal (size (ci), [2, 2]);
%! assert_equal (ci(1,1), 2.486679110991696,  1e-10);
%! assert_equal (ci(1,2), 3.436164115360526, 1e-10);
%! assert_equal (ci(2,1), -1.027166590567854, 1e-10);
%! assert_equal (ci(2,2), -0.967330908318718, 1e-10);
%! assert_equal ((ci(:,1) + ci(:,2)) / 2, m.Coefficients.Estimate, 1e-10);
%! assert_equal (ci(:,2) - ci(:,1), 2 * t * m.Coefficients.SE, 1e-10);

%!test
%! m  = fitlm (X, y, 'interactions');
%! ci = coefCI (m);
%! t  = tinv (0.975, m.DFE);
%! assert_equal (size (ci), [4, 2]);
%! assert_equal (ci(1,1), -0.201030907566802, 1e-10);
%! assert_equal (ci(1,2),  0.516312363642881, 1e-10);
%! assert_equal ((ci(:,1) + ci(:,2)) / 2, m.Coefficients.Estimate, 1e-10);
%! assert_equal (ci(:,2) - ci(:,1), 2 * t * m.Coefficients.SE, 1e-10);

%!test
%! ## constant model (1 coefficient)
%! m  = fitlm (X, y, 'constant');
%! ci = coefCI (m);
%! t  = tinv (0.975, m.DFE);
%! assert_equal (size (ci), [1, 2]);
%! assert_equal ((ci(1,1) + ci(1,2)) / 2, m.Coefficients.Estimate, 1e-10);
%! assert_equal (ci(1,2) - ci(1,1), 2 * t * m.Coefficients.SE, 1e-10);

%!test
%! m  = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! ci = coefCI (m);
%! t  = tinv (0.975, m.DFE);
%! assert_equal (size (ci), [3, 2]);
%! assert_equal (ci(1,1), -0.355978167660141, 1e-10);
%! assert_equal (ci(1,2),  0.516619434992026, 1e-10);
%! assert_equal ((ci(:,1) + ci(:,2)) / 2, m.Coefficients.Estimate, 1e-10);
%! assert_equal (ci(:,2) - ci(:,1), 2 * t * m.Coefficients.SE, 1e-10);

%!test
%! ## rank-deficient: dropped rows give [0,0], active rows are finite
%! m    = fitlm ([ones(n,1), X, X(:,1)+X(:,2)], y);
%! ci   = coefCI (m);
%! drop = find (m.Coefficients.SE == 0);
%! assert_equal (size (ci), [5, 2]);
%! assert_equal (all (all (ci(drop, :) == 0)), true);
%! assert_equal (all (all (isfinite (ci(setdiff (1:5, drop'), :)))), true);

%!test
%! m  = fitlm ([1;1;1;2;2;2;3;3;3], [2.1;2.3;1.9;4.1;3.9;4.2;6.3;5.8;6.1], ...
%!             'linear', 'CategoricalVars', 1);
%! ci = coefCI (m);
%! assert_equal (size (ci), [3, 2]);
%! assert_equal (ci(1,1), 1.80971256321669,  1e-10);
%! assert_equal (ci(1,2), 2.3902874367833,   1e-10);
%! assert_equal (ci(2,1), 1.55613823658119,  1e-10);
%! assert_equal (ci(2,2), 2.37719509675214,  1e-10);
%! assert_equal (ci(3,1), 3.55613823658119,  1e-10);
%! assert_equal (ci(3,2), 4.37719509675214,  1e-10);

%!test
%! [p, F, r] = coefTest (mdl);
%! assert_equal (size (p), [1, 1]);
%! assert_equal (class (p), 'double');
%! assert_equal (p >= 0 && p <= 1, true);
%! assert_equal (F >= 0, true);
%! assert_equal (p, 9.489880832170599e-28, -1e-8);
%! assert_equal (F, 1.283149098426142e+04, -1e-8);
%! assert_equal (r, 2);

%!test
%! ## formula identity
%! [p, F] = coefTest (mdl);
%! k   = mdl.NumCoefficients;
%! H0  = [zeros(k-1, 1), eye(k-1)];
%! b   = mdl.Coefficients.Estimate;
%! V   = mdl.CoefficientCovariance;
%! Hb  = H0 * b;
%! Fm  = (Hb' * ((H0 * V * H0') \ Hb)) / (k - 1);
%! pm  = betainc (mdl.DFE / (mdl.DFE + (k-1) * Fm), mdl.DFE/2, (k-1)/2);
%! assert_equal (F, Fm, -1e-10);
%! assert_equal (p, pm, -1e-10);

%!test
%! ## explicit H matches default
%! k     = mdl.NumCoefficients;
%! H_exp = [zeros(k-1, 1), eye(k-1)];
%! [p1, F1, r1] = coefTest (mdl);
%! [p2, F2, r2] = coefTest (mdl, H_exp);
%! assert_equal (p2, p1, -1e-10);
%! assert_equal (F2, F1, -1e-10);
%! assert_equal (r2, size (H_exp, 1));

%!test
%! ## pinned single and joint H
%! [p1, F1, r1] = coefTest (mdl, [1 0 0]);
%! assert_equal (p1, 0.314859866747774, -1e-8);
%! assert_equal (F1, 1.072634101844537, -1e-8);
%! assert_equal (r1, 1);
%! [p2, F2, r2] = coefTest (mdl, [0 1 0]);
%! assert_equal (p2, 8.937794169018252e-05, -1e-8);
%! assert_equal (F2, 25.985840929474932, -1e-8);
%! assert_equal (r2, 1);
%! [p3, F3, r3] = coefTest (mdl, [0 0 1]);
%! assert_equal (p3, 8.656938305821102e-19, -1e-8);
%! assert_equal (F3, 1.849410599855684e+03, -1e-8);
%! assert_equal (r3, 1);
%! [pm, Fm, rm] = coefTest (mdl, [0 1 0; 0 0 1]);
%! assert_equal (pm, 9.489880832170599e-28, -1e-8);
%! assert_equal (Fm, 1.283149098426142e+04, -1e-8);
%! assert_equal (rm, 2);

%!test
%! ## trivial hypothesis and C=0
%! b        = mdl.Coefficients.Estimate;
%! [p0, F0] = coefTest (mdl, [0 1 0], b(2));
%! assert_equal (F0 < 1e-12, true);
%! assert_equal (p0, 1, 1e-10);
%! [pa, Fa] = coefTest (mdl, [0 1 0], 0);
%! [pb, Fb] = coefTest (mdl, [0 1 0]);
%! assert_equal (pa, pb, -1e-10);
%! assert_equal (Fa, Fb, -1e-10);

%!test
%! ## H with C
%! [pc, Fc, rc] = coefTest (mdl, [0 1 0; 0 0 1], [1.5; -1.0]);
%! assert_equal (pc, 2.833788304242915e-09, -1e-8);
%! assert_equal (Fc, 77.603887650386312, -1e-8);
%! assert_equal (rc, 2);
%! [pr, Fr] = coefTest (mdl, [0 1 0; 0 0 1], [1.5, -1.0]);
%! assert_equal (pr, pc, -1e-10);
%! assert_equal (Fr, Fc, -1e-10);
%! [ps, Fs] = coefTest (mdl, [0 1 0], 1.5);
%! assert_equal (ps, 0.056184159363707, -1e-8);
%! assert_equal (Fs, 4.199865537706047, -1e-8);

%!test
%! ## no-intercept model
%! m = fitlm (X, y, 'Intercept', false);
%! [p, F, r] = coefTest (m);
%! assert_equal (p, 6.060655830723051e-32, -1e-8);
%! assert_equal (F, 2.646694317541346e+04, -1e-8);
%! assert_equal (r, m.NumCoefficients);
%! [p2, F2] = coefTest (m, eye (m.NumCoefficients));
%! assert_equal (p2, p, -1e-10);
%! assert_equal (F2, F, -1e-10);

%!test
%! ## interaction model
%! m = fitlm (X, y, 'interactions');
%! [p, F, r] = coefTest (m);
%! assert_equal (p, 1.164196605688161e-25, -1e-8);
%! assert_equal (F, 8.107508574885546e+03, -1e-8);
%! assert_equal (r, m.NumCoefficients - 1);
%! assert_equal (r != m.NumPredictors, true);

%!test
%! ## weighted model
%! m = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! [p, F, r] = coefTest (m);
%! assert_equal (p, 1.481920976389473e-27, -1e-8);
%! assert_equal (F, 1.217557180481257e+04, -1e-8);
%! assert_equal (r, 2);
%! assert_equal (p, m.ModelFitVsNullModel.Pvalue, -1e-8);

%!test
%! ## categorical model
%! m = fitlm ([1;1;1;2;2;2;3;3;3], [2.1;2.3;1.9;4.1;3.9;4.2;6.3;5.8;6.1], ...
%!            'linear', 'CategoricalVars', 1);
%! [p, F, r] = coefTest (m);
%! assert_equal (p, 1.197590680415813e-06, -1e-8);
%! assert_equal (F, 2.795000000000035e+02, -1e-8);
%! assert_equal (r, 2);
%! [p1, F1] = coefTest (m, [1 0 0]);
%! assert_equal (F1, 3.133421052631613e+02, -1e-8);
%! assert_equal (p1, 2.087464608380450e-06, -1e-8);
%! [p2, F2] = coefTest (m, [0 1 0]);
%! assert_equal (F2, 1.374078947368438e+02, -1e-8);
%! assert_equal (p2, 2.325514143662469e-05, -1e-8);
%! [p3, F3] = coefTest (m, [0 0 1]);
%! assert_equal (F3, 5.589868421052698e+02, -1e-8);
%! assert_equal (p3, 3.757733067786492e-07, -1e-8);

%!test
%! ## constant model
%! m = fitlm (X, y, 'constant');
%! [p, F, r] = coefTest (m);
%! assert_equal (p, 0.000239936408695073, -1e-8);
%! assert_equal (F, 20.3359164947506, -1e-8);
%! assert_equal (r, 1);

%!test
%! ## rank-deficient model
%! m    = fitlm ([ones(n,1), X, X(:,1)+X(:,2)], y);
%! [p, F] = coefTest (m);
%! assert_equal (isnan (p), true);
%! assert_equal (isnan (F), true);
%! drop = find (m.Coefficients.SE == 0);
%! keep = setdiff (2:m.NumCoefficients, drop');
%! H    = zeros (numel (keep), m.NumCoefficients);
%! for i = 1:numel (keep)
%!   H(i, keep(i)) = 1;
%! endfor
%! [p2, F2, r2] = coefTest (m, H);
%! assert_equal (p2, 6.70657058643085e-30, -1e-8);
%! assert_equal (F2, 17716.1864263456, -1e-8);
%! assert_equal (r2, numel (keep));

%!test
%! p = dwtest (mdl);
%! assert_equal (size (p),  [1, 1]);
%! assert_equal (class (p), 'double');
%! [p, DW] = dwtest (mdl);
%! assert_equal (size (DW), [1, 1]);
%! assert_equal (p  >= 0 && p  <= 1, true);
%! assert_equal (DW >= 0 && DW <= 4, true);
%! assert_equal (p,  4.702593821571290e-04, -1e-6);
%! assert_equal (DW, 0.870000704251173, 1e-12);

%!test
%! [p1, DW1] = dwtest (mdl);
%! [p2, DW2] = dwtest (mdl, 'exact', 'both');
%! assert_equal (p1, p2, 1e-14);
%! assert_equal (DW1, DW2, 1e-14);

%!test
%! ## DW is the same for all method and tail options
%! [~, d1] = dwtest (mdl, 'exact', 'both');
%! [~, d2] = dwtest (mdl, 'exact', 'right');
%! [~, d3] = dwtest (mdl, 'exact', 'left');
%! [~, d4] = dwtest (mdl, 'approximate', 'both');
%! [~, d5] = dwtest (mdl, 'approximate', 'right');
%! [~, d6] = dwtest (mdl, 'approximate', 'left');
%! assert_equal (d1, 0.870000704251173, 1e-12);
%! assert_equal (d2, 0.870000704251173, 1e-12);
%! assert_equal (d3, 0.870000704251173, 1e-12);
%! assert_equal (d4, 0.870000704251173, 1e-12);
%! assert_equal (d5, 0.870000704251173, 1e-12);
%! assert_equal (d6, 0.870000704251173, 1e-12);

%!test
%! ## one-sided p-values sum to 1 and two-sided equals twice the smaller
%! pb = dwtest (mdl, 'exact', 'both');
%! pr = dwtest (mdl, 'exact', 'right');
%! pl = dwtest (mdl, 'exact', 'left');
%! assert_equal (pr + pl, 1, 1e-12);
%! assert_equal (pb, 4.702593821571290e-04, 1e-12);

%!test
%! ## all six method and tail combinations pinned
%! assert_equal (dwtest (mdl, 'exact', 'both'), 4.702593821571290e-04, -1e-6);
%! assert_equal (dwtest (mdl, 'exact', 'right'), 2.351296910785645e-04, -1e-6);
%! assert_equal (dwtest (mdl, 'exact', 'left'), 0.999764870308921, -1e-6);
%! assert_equal (dwtest (mdl, 'approximate', 'both'), 0.001058795514879, -1e-6);
%! assert_equal (dwtest (mdl, 'approximate', 'right'), 5.293977574395035e-04, -1e-6);
%! assert_equal (dwtest (mdl, 'approximate', 'left'), 0.999470602242560, -1e-6);

%!test
%! ## no-intercept model
%! m = fitlm (X, y, 'Intercept', false);
%! [p, DW] = dwtest (m, 'exact', 'both');
%! assert_equal (DW, 0.841468411374128, 1e-12);
%! assert_equal (p, 0.001402191159200, -1e-6);
%! assert_equal (dwtest (m, 'exact', 'right'), 7.010955795999754e-04, -1e-6);
%! assert_equal (dwtest (m, 'approximate', 'right'), 0.001350534002321, -1e-6);

%!test
%! ## weighted model
%! m = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! [p, DW] = dwtest (m, 'exact', 'both');
%! assert_equal (DW, 0.871162354803032, 1e-12);
%! assert_equal (p, 4.771641146603785e-04, -1e-6);
%! assert_equal (dwtest (m, 'exact', 'right'), 2.385820573301892e-04, -1e-6);
%! assert_equal (dwtest (m, 'approximate', 'right'), 5.346779629058873e-04, -1e-6);

%!test
%! ## positive autocorrelation model
%! m = fitlm ((1:n)'/n, sin (pi * (1:n)'/n));
%! [~, DW] = dwtest (m, 'exact', 'both');
%! pr = dwtest (m, 'exact', 'right');
%! pl = dwtest (m, 'exact', 'left');
%! assert_equal (DW, 0.118112272685229, 1e-10);
%! assert_equal (DW < 1, true);
%! assert_equal (pr < pl, true);
%! assert_equal (pr < 1e-10, true);

%!test
%! ## negative autocorrelation model
%! m = fitlm ((1:n)'/n, 2*(1:n)'/n + repmat ([1; -1], n/2, 1));
%! [pb, DW] = dwtest (m, 'exact', 'both');
%! pl = dwtest (m, 'exact', 'left');
%! pr = dwtest (m, 'exact', 'right');
%! assert_equal (pb, 4.205713999283489e-09, 1e-10);
%! assert_equal (DW, 3.825974025974026, 1e-10);
%! assert_equal (DW > 2, true);
%! assert_equal (pl < pr, true);
%! assert_equal (pb, 2 * pl, 1e-10);
%! assert_equal (pb < 1e-7, true);

%!test
%! m = addTerms (mdl, 'x1:x2');
%! assert_equal (isa (m, 'LinearModel'), true);
%! assert_equal (mdl.NumCoefficients, 3);
%! assert_equal (m.NumCoefficients, 4);
%! assert_equal (m.NumPredictors, 2);
%! assert_equal (m.NumObservations, 20);
%! assert_equal (m.DFE, 16);
%! assert_equal (m.Coefficients.Estimate(1), 0.157640728038039, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), 2.085426803117909, -1e-8);
%! assert_equal (m.Coefficients.Estimate(3), -0.929682701072813, -1e-8);
%! assert_equal (m.Coefficients.Estimate(4), -0.031208018255475, -1e-8);
%! assert_equal (m.Coefficients.SE(1), 0.169192291625763, -1e-8);
%! assert_equal (m.Coefficients.SE(2), 1.361534257888685, -1e-8);
%! assert_equal (m.Coefficients.SE(3), 0.148744319911833, -1e-8);
%! assert_equal (m.Coefficients.SE(4), 0.0932669056882381, -1e-8);
%! assert_equal (m.Coefficients.tStat(1), 0.931725237144526, -1e-8);
%! assert_equal (m.Coefficients.tStat(2), 1.531674132351069, -1e-8);
%! assert_equal (m.Coefficients.tStat(3), -6.250206405353000, -1e-8);
%! assert_equal (m.Coefficients.tStat(4), -0.334609774230031, -1e-8);
%! assert_equal (m.Coefficients.pValue(1), 0.365325503492671, -1e-8);
%! assert_equal (m.Coefficients.pValue(2), 0.145134783727025, -1e-8);
%! assert_equal (m.Coefficients.pValue(3), 1.159217784590233e-05, -1e-8);
%! assert_equal (m.Coefficients.pValue(4), 0.742265736761240, -1e-8);
%! assert_equal (m.SSE, 0.383859187927621, -1e-8);
%! assert_equal (m.MSE, 0.023991199245515, -1e-8);
%! assert_equal (m.RMSE, 0.154890926930905, -1e-8);
%! assert_equal (m.Rsquared.Ordinary, 0.999342606032059, -1e-8);
%! assert_equal (m.Rsquared.Adjusted, 0.999219344663070, -1e-8);
%! assert_equal (m.LogLikelihood, 11.153346988927943, -1e-8);
%! assert_equal (m.ModelFitVsNullModel.Fstat, 8.107508574898859e+03, -1e-6);
%! assert_equal (m.ModelFitVsNullModel.Pvalue, 1.164196605672873e-25, -1e-6);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.CoefficientNames{3}, 'x2');
%! assert_equal (m.CoefficientNames{4}, 'x1:x2');

%!test
%! ## x1*x2 crossing gives same result as x1:x2 when main effects exist
%! m = addTerms (mdl, 'x1*x2');
%! assert_equal (m.NumCoefficients, 4);
%! assert_equal (m.DFE, 16);
%! assert_equal (m.SSE, 0.383859187927621, -1e-8);
%! assert_equal (m.Coefficients.Estimate(1), 0.157640728038039, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), 2.085426803117909, -1e-8);
%! assert_equal (m.Coefficients.Estimate(3), -0.929682701072813, -1e-8);
%! assert_equal (m.Coefficients.Estimate(4), -0.031208018255475, -1e-8);
%! assert_equal (m.CoefficientNames{4}, 'x1:x2');

%!test
%! m = addTerms (mdl, 'x1 + x1:x2');
%! assert_equal (m.NumCoefficients, 4);
%! assert_equal (m.DFE, 16);
%! assert_equal (m.SSE, 0.383859187927621, -1e-8);
%! assert_equal (m.Coefficients.Estimate(1), 0.157640728038039, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), 2.085426803117909, -1e-8);
%! assert_equal (m.Coefficients.Estimate(3), -0.929682701072813, -1e-8);
%! assert_equal (m.Coefficients.Estimate(4), -0.031208018255475, -1e-8);

%!test
%! ## adding existing term returns equivalent model
%! ws = warning ('off', 'all');
%! m  = addTerms (mdl, 'x1');
%! warning (ws);
%! assert_equal (m.NumCoefficients, 3);
%! assert_equal (m.DFE, 17);
%! assert_equal (m.Coefficients.Estimate(1), 0.116188677790207, 1e-7);
%! assert_equal (m.Coefficients.Estimate(2), 2.508451490570863, 1e-7);
%! assert_equal (m.Coefficients.Estimate(3), -0.978835329825186, 1e-7);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.CoefficientNames{3}, 'x2');

%!test
%! m = addTerms (mdl, 'x2^2');
%! assert_equal (m.NumCoefficients, 4);
%! assert_equal (m.DFE, 16);
%! assert_equal (m.SSE, 0.386103933724971, -1e-8);
%! assert_equal (m.Coefficients.Estimate(1), 0.130152473216993, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), 2.380771990884563, -1e-8);
%! assert_equal (m.Coefficients.Estimate(3), -0.967672823484773, -1e-8);
%! assert_equal (m.Coefficients.Estimate(4), -2.991483322469043e-04, -1e-8);
%! assert_equal (m.Coefficients.SE(1), 0.154974488176692, -1e-8);
%! assert_equal (m.Coefficients.SE(2), 1.071554310049276, -1e-8);
%! assert_equal (m.Coefficients.SE(3), 0.085801310569364, -1e-8);
%! assert_equal (m.Coefficients.SE(4), 0.002211890858232, -1e-8);
%! assert_equal (m.CoefficientNames{4}, 'x2^2');

%!test
%! ## x1^2 rank-deficient: DFE unchanged SE zero for dropped term
%! m = addTerms (mdl, 'x1^2');
%! assert_equal (m.NumCoefficients, 4);
%! assert_equal (m.DFE, 17);
%! assert_equal (m.SSE, 0.386545331386823, -1e-8);
%! assert_equal (m.Coefficients.Estimate(4), 0);
%! assert_equal (m.Coefficients.SE(4), 0);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.CoefficientNames{3}, 'x2');
%! assert_equal (m.CoefficientNames{4}, 'x1^2');


%!test
%! ## numeric matrix [1,1,0] same as string x1:x2
%! m = addTerms (mdl, [1, 1, 0]);
%! assert_equal (m.NumCoefficients, 4);
%! assert_equal (m.DFE, 16);
%! assert_equal (m.SSE, 0.383859187927621, -1e-8);
%! assert_equal (m.Coefficients.Estimate(1), 0.157640728038039, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), 2.085426803117909, -1e-8);
%! assert_equal (m.Coefficients.Estimate(3), -0.929682701072813, -1e-8);
%! assert_equal (m.Coefficients.Estimate(4), -0.031208018255475, -1e-8);

%!test
%! ## numeric matrix [1,1] auto-padded to [1,1,0]
%! m = addTerms (mdl, [1, 1]);
%! assert_equal (m.NumCoefficients, 4);
%! assert_equal (m.DFE, 16);
%! assert_equal (m.SSE, 0.383859187927621, -1e-8);
%! assert_equal (m.Coefficients.Estimate(1), 0.157640728038039, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), 2.085426803117909, -1e-8);
%! assert_equal (m.Coefficients.Estimate(3), -0.929682701072813, -1e-8);
%! assert_equal (m.Coefficients.Estimate(4), -0.031208018255475, -1e-8);

%!test
%! m = addTerms (mdl, [1, 1, 0; 0, 2, 0]);
%! assert_equal (m.NumCoefficients, 5);
%! assert_equal (m.DFE, 15);
%! assert_equal (m.SSE, 0.315784637443501, -1e-8);
%! assert_equal (m.CoefficientNames{4}, 'x1:x2');
%! assert_equal (m.CoefficientNames{5}, 'x2^2');

%!test
%! mc = fitlm (X, y, 'constant');
%! m  = addTerms (mc, 'x1');
%! assert_equal (m.NumCoefficients, 2);
%! assert_equal (m.DFE, 18);
%! assert_equal (m.Coefficients.Estimate(1), 3.884704697617172, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), -18.047090435758047, -1e-8);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');

%!test
%! ## step from constant to full linear model
%! mc  = fitlm (X, y, 'constant');
%! mc1 = addTerms (mc, 'x1');
%! mc2 = addTerms (mc1, 'x2');
%! assert_equal (mc2.NumCoefficients, 3);
%! assert_equal (mc2.DFE, 17);
%! assert_equal (mc2.Coefficients.Estimate(1), 0.116188677790207, 1e-7);
%! assert_equal (mc2.Coefficients.Estimate(2), 2.508451490570863, 1e-7);
%! assert_equal (mc2.Coefficients.Estimate(3), -0.978835329825186, 1e-7);
%! assert_equal (mc2.CoefficientNames{1}, '(Intercept)');
%! assert_equal (mc2.CoefficientNames{2}, 'x1');
%! assert_equal (mc2.CoefficientNames{3}, 'x2');

%!test
%! ## adding intercept to no-intercept model
%! mni = fitlm (X, y, 'Intercept', false);
%! m   = addTerms (mni, '1');
%! assert_equal (m.NumCoefficients, 3);
%! assert_equal (m.DFE, 17);
%! assert_equal (m.Coefficients.Estimate(1), 0.116188677790207, 1e-7);
%! assert_equal (m.Coefficients.Estimate(2), 2.508451490570863, 1e-7);
%! assert_equal (m.Coefficients.Estimate(3), -0.978835329825186, 1e-7);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.CoefficientNames{3}, 'x2');

%!test
%! ## weighted model weights preserved
%! mw = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! m  = addTerms (mw, 'x1:x2');
%! assert_equal (m.NumCoefficients, 4);
%! assert_equal (m.DFE, 16);
%! assert_equal (m.SSE, 0.019230053719402, -1e-8);
%! assert_equal (m.Coefficients.Estimate(1), -0.122645849510537, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), 4.125799311652051, -1e-8);
%! assert_equal (m.Coefficients.Estimate(3), -1.128467507844852, -1e-8);
%! assert_equal (m.Coefficients.Estimate(4), 0.081921546574140, -1e-8);
%! assert_equal (m.Coefficients.SE(1), 0.354926338845420, -1e-8);
%! assert_equal (m.Coefficients.SE(2), 2.205951720932299, -1e-8);
%! assert_equal (m.Coefficients.SE(3), 0.205033663966776, -1e-8);
%! assert_equal (m.Coefficients.SE(4), 0.115523382309399, -1e-8);
%! assert_equal (m.CoefficientNames{4}, 'x1:x2');

%!test
%! ## excluded observations preserved
%! me = fitlm (X, y, 'Exclude', [1, 2]);
%! m  = addTerms (me, 'x1:x2');
%! assert_equal (m.NumObservations, 18);
%! assert_equal (m.DFE, 14);
%! assert_equal (m.NumCoefficients, 4);
%! assert_equal (m.Coefficients.Estimate(1), -0.345521184099998, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), 5.139185607268283, -1e-8);
%! assert_equal (m.Coefficients.Estimate(3), -1.198851436671170, -1e-8);
%! assert_equal (m.Coefficients.Estimate(4), 0.112619530301200, -1e-8);
%! assert_equal (m.CoefficientNames{4}, 'x1:x2');

%!test
%! ## remove two predictors by string from a 4-predictor model
%! Xh = [7 26 6 60; 1 29 15 52; 11 56 8 20; 11 31 8 47; 7 52 6 33; ...
%!        11 55 9 22; 3 71 17 6; 1 31 22 44; 2 54 18 22; 21 47 4 26; ...
%!        1 40 23 34; 11 66 9 12; 10 68 8 12];
%! yh = [78.5;74.3;104.3;87.6;95.9;109.2;102.7;72.5;93.1;115.9;83.8;113.3;109.4];
%! m  = removeTerms (fitlm (Xh, yh), 'x3 + x4');
%! assert_equal (m.NumCoefficients, 3);
%! assert_equal (m.NumEstimatedCoefficients, 3);
%! assert_equal (m.DFE, 10);
%! assert_equal (m.NumObservations, 13);
%! assert_equal (m.NumVariables, 5);
%! assert_equal (m.Coefficients.Estimate(1), 52.577348882089481, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), 1.468305742215555, -1e-8);
%! assert_equal (m.Coefficients.Estimate(3), 0.662250491274645, -1e-8);
%! assert_equal (m.Coefficients.SE(1), 2.286174334503340, -1e-8);
%! assert_equal (m.Coefficients.SE(2), 0.121300923606266, -1e-8);
%! assert_equal (m.Coefficients.SE(3), 0.045854721468522, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.Coefficients.tStat(1), 22.997961305305111, -1e-8);
%! assert_equal (m.Coefficients.tStat(2), 12.104654264476748, -1e-8);
%! assert_equal (m.Coefficients.tStat(3), 14.442362096327519, -1e-8);
%! assert_equal (m.Coefficients.pValue(1), 5.456570901490983e-10, -1e-7);
%! assert_equal (m.Coefficients.pValue(2), 2.692212179685427e-07, -1e-8);
%! assert_equal (m.Coefficients.pValue(3), 5.028960315638413e-08, -1e-8);
%! assert_equal (m.SSE, 57.904483176113658, -1e-8);
%! assert_equal (m.RMSE, 2.40633503852047, -1e-8);
%! assert_equal (m.MSE, 5.790448317611299, 1e-12);
%! assert_equal (m.SST, 2.715763076923078e+03, 1e-8);
%! assert_equal (m.Rsquared.Ordinary, 0.978678374535632, -1e-8);
%! assert_equal (m.Rsquared.Adjusted, 0.974414049442758, -1e-8);
%! assert_equal (size (m.CoefficientCovariance), [3, 3]);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.CoefficientNames{3}, 'x2');
%! assert_equal (m.Formula.HasIntercept, true);
%! assert_equal (m.Formula.LinearPredictor, '1 + x1 + x2');
%! assert_equal (height (m.Diagnostics), 13);
%! assert_equal (sum (m.Diagnostics.Leverage), 3, 1e-10);
%! assert_equal (m.Residuals.Raw, yh - m.Fitted, 1e-10);

%!test
%! m = removeTerms (mdl, 'x2');
%! assert_equal (m.NumCoefficients, 2);
%! assert_equal (m.NumEstimatedCoefficients, 2);
%! assert_equal (m.DFE, 18);
%! assert_equal (m.NumObservations, 20);
%! assert_equal (m.Coefficients.Estimate(1), 3.88470469761717, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), -18.047090435758, -1e-8);
%! assert_equal (m.Coefficients.SE(2), 1.19086428900602, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.Coefficients.tStat(2), -15.1546155194741, -1e-8);
%! assert_equal (m.SSE, 42.4383708132815, -1e-8);
%! assert_equal (m.SST, 583.910420002346, -1e-8);
%! assert_equal (m.Rsquared.Ordinary, 0.927320408474452, -1e-8);
%! assert_equal (size (m.CoefficientCovariance), [2, 2]);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.Formula.HasIntercept, true);
%! assert_equal (m.Formula.LinearPredictor, '1 + x1');
%! assert_equal (height (m.Diagnostics), 20);
%! assert_equal (m.Residuals.Raw, y - m.Fitted, 1e-10);

%!test
%! ## removing the intercept via string '1'
%! m = removeTerms (mdl, '1');
%! assert_equal (m.NumCoefficients, 2);
%! assert_equal (m.NumEstimatedCoefficients, 2);
%! assert_equal (m.DFE, 18);
%! assert_equal (m.NumObservations, 20);
%! assert_equal (m.Formula.HasIntercept, false);
%! assert_equal (m.Formula.LinearPredictor, 'x1 + x2');
%! assert_equal (m.Coefficients.Estimate(1), 2.96142161317611, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), -0.997248749443286, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.SSE, 0.410934843407688, -1e-8);
%! assert_equal (size (m.CoefficientCovariance), [2, 2]);
%! assert_equal (m.CoefficientNames{1}, 'x1');
%! assert_equal (m.CoefficientNames{2}, 'x2');
%! assert_equal (! any (strcmp (m.CoefficientNames, '(Intercept)')), true);
%! assert_equal (height (m.Diagnostics), 20);

%!test
%! ## removing both predictors leaves only the intercept
%! m = removeTerms (mdl, 'x1 + x2');
%! assert_equal (m.NumCoefficients, 1);
%! assert_equal (m.NumEstimatedCoefficients, 1);
%! assert_equal (m.DFE, 19);
%! assert_equal (m.NumObservations, 20);
%! assert_equal (m.Formula.HasIntercept, true);
%! assert_equal (m.Formula.LinearPredictor, '1');
%! assert_equal (m.Coefficients.Estimate(1), -5.5900177811558, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.SSE, 583.910420002346, -1e-8);
%! assert_equal (m.SST, 583.910420002346, -1e-8);
%! assert_equal (m.SSR, 0, 1e-20);
%! assert_equal (size (m.CoefficientCovariance), [1, 1]);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (height (m.Diagnostics), 20);
%! assert_equal (sum (m.Diagnostics.Leverage), 1, 1e-10);
%! assert_equal (m.Residuals.Raw, y - m.Fitted, 1e-10);

%!test
%! ws = warning ('off', 'all');
%! m = removeTerms (mdl, 'x1:x2');
%! warning (ws);
%! assert_equal (m.NumCoefficients, mdl.NumCoefficients);
%! assert_equal (m.NumEstimatedCoefficients, mdl.NumEstimatedCoefficients);
%! assert_equal (m.DFE, mdl.DFE);
%! assert_equal (m.SSE, mdl.SSE, 1e-15);
%! assert_equal (m.SSR, mdl.SSR, 1e-15);
%! assert_equal (m.SST, mdl.SST, 1e-15);
%! assert_equal (m.RMSE, mdl.RMSE, 1e-15);
%! assert_equal (m.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-15);
%! assert_equal (m.Coefficients.SE, mdl.Coefficients.SE, 1e-15);
%! assert_equal (m.CoefficientCovariance, mdl.CoefficientCovariance, 1e-15);
%! assert_equal (isequal (m.CoefficientNames, mdl.CoefficientNames), true);
%! assert_equal (m.Formula.LinearPredictor, mdl.Formula.LinearPredictor);
%! assert_equal (m.Formula.HasIntercept, mdl.Formula.HasIntercept);

%!test
%! m = removeTerms (mdl, [0 1 0]);
%! assert_equal (m.NumCoefficients, 2);
%! assert_equal (m.NumEstimatedCoefficients, 2);
%! assert_equal (m.DFE, 18);
%! assert_equal (m.NumObservations, 20);
%! assert_equal (m.Coefficients.Estimate(1), 3.88470469761717, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), -18.047090435758, -1e-8);
%! assert_equal (m.Coefficients.SE(2), 1.19086428900602, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.Coefficients.tStat(2), -15.1546155194741, -1e-8);
%! assert_equal (m.SSE, 42.4383708132815, -1e-8);
%! assert_equal (m.SST, 583.910420002346, -1e-8);
%! assert_equal (m.Rsquared.Ordinary, 0.927320408474452, -1e-8);
%! assert_equal (size (m.CoefficientCovariance), [2, 2]);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.Formula.HasIntercept, true);
%! assert_equal (m.Formula.LinearPredictor, '1 + x1');
%! assert_equal (height (m.Diagnostics), 20);
%! assert_equal (sum (m.Diagnostics.Leverage), 2, 1e-10);

%!test
%! ## auto-padded matrix [0 1] gives identical result to [0 1 0]
%! m = removeTerms (mdl, [0 1]);
%! assert_equal (m.NumCoefficients, 2);
%! assert_equal (m.NumEstimatedCoefficients, 2);
%! assert_equal (m.DFE, 18);
%! assert_equal (m.Coefficients.Estimate(1), 3.88470469761717, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), -18.047090435758, -1e-8);
%! assert_equal (m.Coefficients.SE(2), 1.19086428900602, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.SSE, 42.4383708132815, -1e-8);
%! assert_equal (m.SST, 583.910420002346, -1e-8);
%! assert_equal (m.Rsquared.Ordinary, 0.927320408474452, -1e-8);
%! assert_equal (size (m.CoefficientCovariance), [2, 2]);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.Formula.LinearPredictor, '1 + x1');
%! assert_equal (height (m.Diagnostics), 20);
%! assert_equal (sum (m.Diagnostics.Leverage), 2, 1e-10);

%!test
%! ## multi-row matrix removes two terms same as the string form
%! Xh = [7 26 6 60; 1 29 15 52; 11 56 8 20; 11 31 8 47; 7 52 6 33; ...
%!        11 55 9 22; 3 71 17 6; 1 31 22 44; 2 54 18 22; 21 47 4 26; ...
%!        1 40 23 34; 11 66 9 12; 10 68 8 12];
%! yh = [78.5;74.3;104.3;87.6;95.9;109.2;102.7;72.5;93.1;115.9;83.8;113.3;109.4];
%! m = removeTerms (fitlm (Xh, yh), [0 0 1 0 0; 0 0 0 1 0]);
%! assert_equal (m.NumCoefficients, 3);
%! assert_equal (m.NumEstimatedCoefficients, 3);
%! assert_equal (m.DFE, 10);
%! assert_equal (m.NumObservations, 13);
%! assert_equal (m.Coefficients.Estimate(1), 52.577348882089481, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), 1.468305742215555, -1e-8);
%! assert_equal (m.Coefficients.Estimate(3), 0.662250491274645, -1e-8);
%! assert_equal (m.Coefficients.SE(1), 2.286174334503340, -1e-8);
%! assert_equal (m.Coefficients.SE(2), 0.121300923606266, -1e-8);
%! assert_equal (m.Coefficients.SE(3), 0.045854721468522, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.Coefficients.tStat(2), 12.104654264476748, -1e-8);
%! assert_equal (m.Coefficients.tStat(3), 14.442362096327519, -1e-8);
%! assert_equal (m.Coefficients.pValue(2), 2.692212179685427e-07, -1e-8);
%! assert_equal (m.SSE, 57.904483176113658, -1e-8);
%! assert_equal (m.Rsquared.Ordinary, 0.978678374535632, -1e-8);
%! assert_equal (m.Rsquared.Adjusted, 0.974414049442758, -1e-8);
%! assert_equal (size (m.CoefficientCovariance), [3, 3]);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.CoefficientNames{3}, 'x2');
%! assert_equal (m.Formula.HasIntercept, true);
%! assert_equal (m.Formula.LinearPredictor, '1 + x1 + x2');
%! assert_equal (height (m.Diagnostics), 13);
%! assert_equal (sum (m.Diagnostics.Leverage), 3, 1e-10);
%! assert_equal (m.Residuals.Raw, yh - m.Fitted, 1e-10);

%!test
%! ## observation weights carry through to the refitted model
%! w = (1:n)' / sum (1:n);
%! mw = fitlm (X, y, 'Weights', w);
%! m = removeTerms (mw, 'x2');
%! assert_equal (m.NumCoefficients, 2);
%! assert_equal (m.NumEstimatedCoefficients, 2);
%! assert_equal (m.DFE, 18);
%! assert_equal (m.NumObservations, 20);
%! assert_equal (m.Coefficients.Estimate(1), 6.29263960898714, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), -21.5708976231287, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.SSE, 1.41763159723151, -1e-8);
%! assert_equal (size (m.CoefficientCovariance), [2, 2]);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.Formula.HasIntercept, true);
%! assert_equal (m.Formula.LinearPredictor, '1 + x1');
%! assert_equal (m.ObservationInfo.Weights, w, 1e-15);
%! assert_equal (sum (m.ObservationInfo.Weights), 1, 1e-12);
%! assert_equal (height (m.Diagnostics), 20);
%! assert_equal (sum (m.Diagnostics.Leverage), 2, 1e-10);
%! assert_equal (m.SSE != removeTerms (mdl, 'x2').SSE, true);

%!test
%! ## excluded rows are preserved and reduce effective sample size
%! me = fitlm (X, y, 'Exclude', [1, 3]);
%! m = removeTerms (me, 'x2');
%! assert_equal (m.NumObservations, 18);
%! assert_equal (m.DFE, 16);
%! assert_equal (m.NumCoefficients, 2);
%! assert_equal (m.NumEstimatedCoefficients, 2);
%! assert_equal (m.Coefficients.Estimate(1), 4.96609542902066, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), -19.5618050042778, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (size (m.CoefficientCovariance), [2, 2]);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.Formula.LinearPredictor, '1 + x1');
%! assert_equal (m.ObservationInfo.Excluded(1), true);
%! assert_equal (m.ObservationInfo.Excluded(3), true);
%! assert_equal (m.ObservationInfo.Excluded(2), false);
%! assert_equal (m.ObservationInfo.Missing(1), false);
%! assert_equal (height (m.Diagnostics), 20);
%! assert_equal (isnan (m.Fitted(1)), true);
%! assert_equal (isnan (m.Fitted(3)), true);
%! assert_equal (isfinite (m.Fitted(2)), true);

%!test
%! ## removing x2 from a no-intercept model gives one slope term
%! mni = fitlm (X, y, 'Intercept', false);
%! m = removeTerms (mni, 'x2');
%! assert_equal (m.NumCoefficients, 1);
%! assert_equal (m.NumEstimatedCoefficients, 1);
%! assert_equal (m.DFE, 19);
%! assert_equal (m.NumObservations, 20);
%! assert_equal (m.Formula.HasIntercept, false);
%! assert_equal (m.Formula.LinearPredictor, 'x1');
%! assert_equal (m.Coefficients.Estimate(1), -12.362156731928, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.SSE, 112.371951585499, -1e-8);
%! assert_equal (size (m.CoefficientCovariance), [1, 1]);
%! assert_equal (m.CoefficientCovariance(1,1) > 0, true);
%! assert_equal (m.CoefficientNames{1}, 'x1');
%! assert_equal (! any (strcmp (m.CoefficientNames, '(Intercept)')), true);
%! assert_equal (height (m.Diagnostics), 20);
%! assert_equal (all (isfinite (m.Fitted)), true);
%! assert_equal (sum (m.Diagnostics.Leverage), 1, 1e-10);

%!test
%! ## removing the interaction term recovers the plain linear model
%! mi = fitlm (X, y, 'interactions');
%! m = removeTerms (mi, 'x1:x2');
%! assert_equal (m.NumCoefficients, 3);
%! assert_equal (m.NumEstimatedCoefficients, 3);
%! assert_equal (m.DFE, 17);
%! assert_equal (m.NumObservations, 20);
%! assert_equal (m.Coefficients.Estimate(1), 0.116188677790207, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), 2.508451490570863, -1e-8);
%! assert_equal (m.Coefficients.Estimate(3), -0.978835329825186, -1e-8);
%! assert_equal (m.Coefficients.SE(1), 0.112185831, -1e-7);
%! assert_equal (m.Coefficients.SE(2), 0.4920818186, -1e-7);
%! assert_equal (m.Coefficients.SE(3), 0.02276108523, -1e-7);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.Coefficients.tStat(1), 1.035680502, -1e-6);
%! assert_equal (m.Coefficients.tStat(2), 5.097630913, -1e-6);
%! assert_equal (m.Coefficients.tStat(3), -43.00477415, -1e-6);
%! assert_equal (m.SSE, 0.386545331386823, -1e-8);
%! assert_equal (size (m.CoefficientCovariance), [3, 3]);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.CoefficientNames{3}, 'x2');
%! assert_equal (m.Formula.HasIntercept, true);
%! assert_equal (m.Formula.LinearPredictor, '1 + x1 + x2');
%! assert_equal (height (m.Diagnostics), 20);
%! assert_equal (sum (m.Diagnostics.Leverage), 3, 1e-10);

%!test
%! ## removing a quadratic term refits on the remaining terms
%! mq = fitlm (X, y, 'quadratic');
%! m = removeTerms (mq, 'x2^2');
%! assert_equal (m.NumEstimatedCoefficients, 4);
%! assert_equal (m.DFE, 16);
%! assert_equal (m.SSE, 0.383859187927621, -1e-8);
%! assert_equal (size (m.CoefficientCovariance, 1), m.NumCoefficients);
%! assert_equal (size (m.CoefficientCovariance, 2), m.NumCoefficients);
%! assert_equal (m.Formula.HasIntercept, true);
%! assert_equal (height (m.Diagnostics), 20);
%! assert_equal (sum (m.Diagnostics.Leverage), 4, 1e-10);
%! assert_equal (m.SSE >= mq.SSE, true);
%! assert_equal (! any (strcmp (m.CoefficientNames, 'x2^2')), true);

%!test
%! ## star notation removes main effects and interaction in one call
%! mi = fitlm (X, y, 'interactions');
%! m = removeTerms (mi, 'x1*x2');
%! assert_equal (m.NumCoefficients, 1);
%! assert_equal (m.NumEstimatedCoefficients, 1);
%! assert_equal (m.DFE, 19);
%! assert_equal (m.NumObservations, 20);
%! assert_equal (m.Formula.HasIntercept, true);
%! assert_equal (m.Formula.LinearPredictor, '1');
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.Coefficients.Estimate(1), -5.5900177811558, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.SSE, 583.910420002346, -1e-8);
%! assert_equal (m.SST, 583.910420002346, -1e-8);
%! assert_equal (m.SSR, 0, 1e-20);
%! assert_equal (size (m.CoefficientCovariance), [1, 1]);
%! assert_equal (height (m.Diagnostics), 20);
%! assert_equal (sum (m.Diagnostics.Leverage), 1, 1e-10);

%!test
%! ## 3-predictor model: removing one term matches a direct two-predictor fit
%! X3 = [X, sin((1:n)' * pi / n)];
%! y3 = X3 * [3; -1; 2] + 0.1 * cos ((1:n)' * pi / 7);
%! m = removeTerms (fitlm (X3, y3), 'x3');
%! r = fitlm (X, y3);
%! assert_equal (m.NumCoefficients, 3);
%! assert_equal (m.NumEstimatedCoefficients, 3);
%! assert_equal (m.DFE, 17);
%! assert_equal (m.NumObservations, 20);
%! assert_equal (m.Coefficients.Estimate, r.Coefficients.Estimate, 1e-10);
%! assert_equal (m.Coefficients.SE, r.Coefficients.SE, 1e-10);
%! assert_equal (m.Coefficients.tStat, r.Coefficients.tStat, 1e-10);
%! assert_equal (m.Coefficients.pValue, r.Coefficients.pValue, 1e-10);
%! assert_equal (m.SSE, r.SSE, 1e-12);
%! assert_equal (m.SSR, r.SSR, 1e-12);
%! assert_equal (m.SST, r.SST, 1e-12);
%! assert_equal (m.MSE, r.MSE, 1e-12);
%! assert_equal (m.RMSE, r.RMSE, 1e-12);
%! assert_equal (m.Rsquared.Ordinary, r.Rsquared.Ordinary, 1e-12);
%! assert_equal (m.Rsquared.Adjusted, r.Rsquared.Adjusted, 1e-12);
%! assert_equal (m.CoefficientCovariance, r.CoefficientCovariance, 1e-12);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.CoefficientNames{3}, 'x2');
%! assert_equal (m.Formula.HasIntercept, true);
%! assert_equal (m.Formula.LinearPredictor, '1 + x1 + x2');
%! assert_equal (height (m.Diagnostics), 20);
%! assert_equal (sum (m.Diagnostics.Leverage), 3, 1e-10);

%!test
%! ## 3-predictor model: removing two terms matches a direct one-predictor fit
%! X3 = [X, sin((1:n)' * pi / n)];
%! y3 = X3 * [3; -1; 2] + 0.1 * cos ((1:n)' * pi / 7);
%! m = removeTerms (fitlm (X3, y3), 'x2 + x3');
%! r = fitlm (X(:,1), y3);
%! assert_equal (m.NumCoefficients, 2);
%! assert_equal (m.NumEstimatedCoefficients, 2);
%! assert_equal (m.DFE, 18);
%! assert_equal (m.NumObservations, 20);
%! assert_equal (m.Coefficients.Estimate, r.Coefficients.Estimate, 1e-10);
%! assert_equal (m.Coefficients.SE, r.Coefficients.SE, 1e-10);
%! assert_equal (m.Coefficients.tStat, r.Coefficients.tStat, 1e-10);
%! assert_equal (m.Coefficients.pValue, r.Coefficients.pValue, 1e-10);
%! assert_equal (m.SSE, r.SSE, 1e-12);
%! assert_equal (m.SST, r.SST, 1e-12);
%! assert_equal (m.MSE, r.MSE, 1e-12);
%! assert_equal (m.RMSE, r.RMSE, 1e-12);
%! assert_equal (m.Rsquared.Ordinary, r.Rsquared.Ordinary, 1e-12);
%! assert_equal (m.Rsquared.Adjusted, r.Rsquared.Adjusted, 1e-12);
%! assert_equal (m.CoefficientCovariance, r.CoefficientCovariance, 1e-12);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.Formula.HasIntercept, true);
%! assert_equal (m.Formula.LinearPredictor, '1 + x1');
%! assert_equal (height (m.Diagnostics), 20);
%! assert_equal (sum (m.Diagnostics.Leverage), 2, 1e-10);

%!test
%! Xh = [7 26 6 60; 1 29 15 52; 11 56 8 20; 11 31 8 47; 7 52 6 33; ...
%!        11 55 9 22; 3 71 17 6; 1 31 22 44; 2 54 18 22; 21 47 4 26; ...
%!        1 40 23 34; 11 66 9 12; 10 68 8 12];
%! yh = [78.5;74.3;104.3;87.6;95.9;109.2;102.7;72.5;93.1;115.9;83.8;113.3;109.4];
%! m = removeTerms (fitlm (Xh, yh), 'x3 + x4');
%! r = removeTerms (removeTerms (fitlm (Xh, yh), 'x4'), 'x3');
%! assert_equal (r.NumCoefficients, 3);
%! assert_equal (r.DFE, 10);
%! assert_equal (r.SSE, m.SSE, 1e-12);
%! assert_equal (r.SSR, m.SSR, 1e-12);
%! assert_equal (r.SST, m.SST, 1e-12);
%! assert_equal (r.RMSE, m.RMSE, 1e-12);
%! assert_equal (r.Coefficients.Estimate, m.Coefficients.Estimate, 1e-10);
%! assert_equal (r.Coefficients.SE, m.Coefficients.SE, 1e-10);
%! assert_equal (r.Coefficients.tStat, m.Coefficients.tStat, 1e-10);
%! assert_equal (r.Coefficients.pValue, m.Coefficients.pValue, 1e-10);
%! assert_equal (r.CoefficientCovariance, m.CoefficientCovariance, 1e-12);
%! assert_equal (isequal (r.CoefficientNames, m.CoefficientNames), true);
%! assert_equal (r.Rsquared.Ordinary, m.Rsquared.Ordinary, 1e-12);
%! assert_equal (r.Rsquared.Adjusted, m.Rsquared.Adjusted, 1e-12);
%! assert_equal (r.Formula.LinearPredictor, m.Formula.LinearPredictor);

%!test
%! m = removeTerms (mdl, [0 0 0]);
%! r = removeTerms (mdl, '1');
%! assert_equal (m.NumCoefficients, 2);
%! assert_equal (m.Formula.HasIntercept, false);
%! assert_equal (m.Formula.LinearPredictor, 'x1 + x2');
%! assert_equal (m.Coefficients.Estimate(1), 2.96142161317611, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), -0.997248749443286, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.SSE, r.SSE, 1e-15);
%! assert_equal (m.SSR, r.SSR, 1e-15);
%! assert_equal (m.SST, r.SST, 1e-15);
%! assert_equal (m.Coefficients.Estimate, r.Coefficients.Estimate, 1e-15);
%! assert_equal (m.Coefficients.SE, r.Coefficients.SE, 1e-15);
%! assert_equal (m.CoefficientCovariance, r.CoefficientCovariance, 1e-15);
%! assert_equal (isequal (m.CoefficientNames, r.CoefficientNames), true);
%! assert_equal (height (m.Diagnostics), 20);
%! assert_equal (sum (m.Diagnostics.Leverage), 2, 1e-10);

%!test
%! ## removing a categorical predictor drops all its indicator variables at once
%! Xc = [1;1;1;2;2;2;3;3;3];
%! yc = [2.1;2.3;1.9; 4.1;3.9;4.2; 6.3;5.8;6.1];
%! mc = fitlm (Xc, yc, 'linear', 'CategoricalVars', 1);
%! m = removeTerms (mc, 'x1');
%! assert_equal (m.NumCoefficients, 1);
%! assert_equal (m.NumEstimatedCoefficients, 1);
%! assert_equal (m.DFE, 8);
%! assert_equal (m.NumObservations, 9);
%! assert_equal (m.Formula.HasIntercept, true);
%! assert_equal (m.Formula.LinearPredictor, '1');
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.Coefficients.Estimate(1), 4.07777777777778, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.SSR, 0, 1e-20);
%! assert_equal (size (m.CoefficientCovariance), [1, 1]);
%! assert_equal (height (m.Diagnostics), 9);
%! assert_equal (all (isfinite (m.Fitted)), true);
%! assert_equal (sum (m.Diagnostics.Leverage), 1, 1e-10);
%! assert_equal (m.Residuals.Raw, yc - m.Fitted, 1e-10);

%!test
%! ## matrix row removes x4 from hald leaving intercept plus x1 x2 x3
%! Xh = [7 26 6 60; 1 29 15 52; 11 56 8 20; 11 31 8 47; 7 52 6 33; ...
%!        11 55 9 22; 3 71 17 6; 1 31 22 44; 2 54 18 22; 21 47 4 26; ...
%!        1 40 23 34; 11 66 9 12; 10 68 8 12];
%! yh = [78.5;74.3;104.3;87.6;95.9;109.2;102.7;72.5;93.1;115.9;83.8;113.3;109.4];
%! m = removeTerms (fitlm (Xh, yh), [0 0 0 1 0]);
%! assert_equal (m.NumCoefficients, 4);
%! assert_equal (m.NumEstimatedCoefficients, 4);
%! assert_equal (m.DFE, 9);
%! assert_equal (m.NumObservations, 13);
%! assert_equal (m.Coefficients.Estimate(1), 48.1936343180437, -1e-8);
%! assert_equal (m.Coefficients.Estimate(2), 1.69589016748479, -1e-8);
%! assert_equal (m.Coefficients.Estimate(3), 0.656914878270554, -1e-8);
%! assert_equal (m.Coefficients.Estimate(4), 0.250017606680009, -1e-8);
%! assert_equal (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert_equal (m.SSE, 48.1106140726532, -1e-8);
%! assert_equal (size (m.CoefficientCovariance), [4, 4]);
%! assert_equal (m.CoefficientNames{1}, '(Intercept)');
%! assert_equal (m.CoefficientNames{2}, 'x1');
%! assert_equal (m.CoefficientNames{3}, 'x2');
%! assert_equal (m.CoefficientNames{4}, 'x3');
%! assert_equal (m.Formula.HasIntercept, true);
%! assert_equal (m.Formula.LinearPredictor, '1 + x1 + x2 + x3');
%! assert_equal (height (m.Diagnostics), 13);
%! assert_equal (all (isfinite (m.Fitted)), true);
%! assert_equal (sum (m.Diagnostics.Leverage), 4, 1e-10);
%! assert_equal (m.Residuals.Raw, yh - m.Fitted, 1e-10);

%!test
%! ## default call creates a histogram with correct bin count and density
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl);
%! xd = get (h(1), 'XData');
%! yd = get (h(1), 'YData');
%! r = mdl.Residuals.Raw(! isnan (mdl.Residuals.Raw));
%! bw = xd(3,1) - xd(1,1);
%! assert_equal (numel (h), 1);
%! assert_equal (get (h(1), 'type'), 'patch');
%! assert_equal (size (xd, 2) > 0, true);
%! assert_equal (sum (yd(2,:)) * bw, 1, 1e-10);
%! assert_equal (all (yd(1,:) == 0) && all (yd(4,:) == 0), true);
%! assert_equal (get (get (ax, 'xlabel'), 'string'), 'Residuals');
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Probability density');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Histogram of residuals');
%! close (fig);

%!test
%! ## histogram bar color changes when FaceColor is passed
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'histogram', 'FaceColor', [0 1 0]);
%! assert_equal (get (h(1), 'FaceColor'), [0 1 0], 1e-10);
%! close (fig);

%!test
%! ## fitted plot shows residuals against fitted values with a zero reference line
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'fitted');
%! assert_equal (numel (h), 2);
%! assert_equal (get (h(1), 'XData'), mdl.Fitted', 1e-15);
%! assert_equal (get (h(1), 'YData'), mdl.Residuals.Raw', 1e-15);
%! assert_equal (get (h(1), 'LineStyle'), 'none');
%! assert_equal (get (h(1), 'Marker'), 'x');
%! assert_equal (get (h(2), 'YData'), [0 0]);
%! assert_equal (get (h(2), 'LineStyle'), ':');
%! assert_equal (get (get (ax, 'xlabel'), 'string'), 'Fitted values');
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Residuals');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Plot of residuals vs. fitted values');
%! close (fig);

%!test
%! ## custom color applies to data points but leaves the reference line unchanged
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'fitted', 'Color', [1 0 0]);
%! assert_equal (get (h(1), 'Color'), [1 0 0], 1e-10);
%! assert_equal (get (h(2), 'Color'), [0.8510 0.8510 0.8510], 1e-4);
%! assert_equal (get (h(2), 'LineStyle'), ':');
%! close (fig);

%!test
%! ## excluded rows appear as gaps in the fitted plot
%! me = fitlm (X, y, 'Exclude', [3, 8]);
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, me, 'fitted');
%! yd = get (h(1), 'YData');
%! assert_equal (numel (yd), 20);
%! assert_equal (isnan (yd(3)), true);
%! assert_equal (isnan (yd(8)), true);
%! assert_equal (! isnan (yd(1)), true);
%! close (fig);

%!test
%! ## case order plot covers all rows and shows gaps where rows were excluded
%! me = fitlm (X, y, 'Exclude', [2, 5]);
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, me, 'caseorder');
%! xd = get (h(1), 'XData');
%! yd = get (h(1), 'YData');
%! assert_equal (xd, 1:20);
%! assert_equal (isnan (yd(2)), true);
%! assert_equal (isnan (yd(5)), true);
%! assert_equal (! isnan (yd(1)), true);
%! assert_equal (get (h(2), 'YData'), [0 0]);
%! assert_equal (get (get (ax, 'xlabel'), 'string'), 'Row number');
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Residuals');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Case order plot of residuals');
%! close (fig);

%!test
%! ## lagged plot shows each residual against the previous one with two reference lines
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'lagged');
%! r = mdl.Residuals.Raw;
%! assert_equal (numel (h), 3);
%! assert_equal (get (h(1), 'XData'), r(1:end-1)', 1e-15);
%! assert_equal (get (h(1), 'YData'), r(2:end)', 1e-15);
%! assert_equal (get (h(2), 'YData'), [0 0]);
%! assert_equal (get (h(3), 'XData'), [0 0]);
%! assert_equal (get (h(2), 'LineStyle'), ':');
%! assert_equal (get (get (ax, 'xlabel'), 'string'), 'Residual(t-1)');
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Residual(t)');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Plot of residuals vs. lagged residuals');
%! close (fig);

%!test
%! ## probability plot uses sorted active residuals as its data
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'probability');
%! r_s = sort (mdl.Residuals.Raw(! isnan (mdl.Residuals.Raw)));
%! assert_equal (get (h(1), 'XData'), r_s', 1e-15);
%! assert_equal (get (get (ax, 'xlabel'), 'string'), 'Residuals');
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Probability');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Normal probability plot of residuals');
%! close (fig);

%!test
%! ## observed plot connects each point to the reference line with a vertical segment
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'observed');
%! obs = mdl.Variables{:, mdl.ResponseName};
%! assert_equal (numel (h), 3);
%! assert_equal (get (h(1), 'XData'), mdl.Fitted', 1e-15);
%! assert_equal (get (h(1), 'YData'), obs', 1e-15);
%! assert_equal (isequal (get (h(2), 'XData'), get (h(2), 'YData')), true);
%! xd3 = get (h(3), 'XData');
%! assert_equal (numel (xd3), 3 * mdl.NumObservations);
%! assert_equal (sum (isnan (xd3)), mdl.NumObservations);
%! assert_equal (get (get (ax, 'xlabel'), 'string'), 'Fitted values');
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Observed response values');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Plot of observed vs. fitted values');
%! close (fig);

%!test
%! ## symmetry plot measures distance from median in both tails
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'symmetry');
%! r_s = sort (mdl.Residuals.Raw(! isnan (mdl.Residuals.Raw)));
%! med = median (r_s);
%! m = floor (numel (r_s) / 2);
%! x_sym = sort (med - r_s(1:m));
%! y_sym = sort (r_s(end-m+1:end) - med);
%! assert_equal (numel (h), 2);
%! assert_equal (get (h(1), 'XData'), x_sym', 1e-15);
%! assert_equal (get (h(1), 'YData'), y_sym', 1e-15);
%! assert_equal (isequal (get (h(2), 'XData'), get (h(2), 'YData')), true);
%! assert_equal (get (h(2), 'LineStyle'), ':');
%! assert_equal (get (get (ax, 'xlabel'), 'string'), 'Lower tail');
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Upper tail');
%! close (fig);

%!test
%! ## switching to pearson residuals changes plotted values but not x positions
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'fitted', 'ResidualType', 'pearson');
%! assert_equal (get (h(1), 'YData'), mdl.Residuals.Pearson', 1e-15);
%! assert_equal (get (h(1), 'XData'), mdl.Fitted', 1e-15);
%! assert_equal (! isequal (get (h(1), 'YData'), mdl.Residuals.Raw'), true);
%! close (fig);

%!test
%! ## standardized and studentized residuals produce different values
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h1 = plotResiduals (ax, mdl, 'caseorder', 'ResidualType', 'standardized');
%! h2 = plotResiduals (ax, mdl, 'caseorder', 'ResidualType', 'studentized');
%! assert_equal (get (h1(1), 'YData'), mdl.Residuals.Standardized', 1e-15);
%! assert_equal (get (h2(1), 'YData'), mdl.Residuals.Studentized', 1e-15);
%! assert_equal (! isequal (get (h1(1), 'YData'), get (h2(1), 'YData')), true);
%! close (fig);

%!test
%! ## marker style and size apply to data points but not to the reference line
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'fitted', 'Marker', 's', 'MarkerSize', 10);
%! assert_equal (get (h(1), 'Marker'), 's');
%! assert_equal (get (h(1), 'MarkerSize'), 10);
%! assert_equal (get (h(2), 'Marker'), 'none');
%! close (fig);

%!test
%! ## weighted model residuals differ from unweighted residuals in the fitted plot
%! mw = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mw, 'fitted');
%! assert_equal (get (h(1), 'YData'), mw.Residuals.Raw', 1e-15);
%! assert_equal (! isequal (get (h(1), 'YData'), mdl.Residuals.Raw'), true);
%! close (fig);

%!test
%! ## calling without an axes handle plots into the current axes
%! fig = figure ('visible', 'off');
%! h = plotResiduals (mdl, 'fitted');
%! assert_equal (isgraphics (get (h(1), 'Parent'), 'axes'), true);
%! assert_equal (isequal (get (h(1), 'Parent'), gca ()), true);
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl);
%! yd = get (h(1), 'YData');
%! assert_equal (numel (h), 2);
%! assert_equal (get (h(1), 'XData'), 1:n);
%! assert_equal (yd(1), 0.370779220779221, -1e-10);
%! assert_equal (yd(2), 0.245283663704716, -1e-10);
%! assert_equal (yd(3), 0.164718614718615, -1e-10);
%! assert_equal (yd(4), 0.118147641831852, -1e-10);
%! assert_equal (yd(5), 0.0960013670539986, -1e-10);
%! assert_equal (yd(6), 0.0900774663932558, -1e-10);
%! assert_equal (yd(7), 0.0935406698564593, -1e-10);
%! assert_equal (yd(8), 0.100922761449077, -1e-10);
%! assert_equal (yd(9), 0.108122579175211, -1e-10);
%! assert_equal (yd(10), 0.112406015037594, -1e-10);
%! assert_equal (yd(11), 0.112406015037594, -1e-10);
%! assert_equal (yd(12), 0.108122579175211, -1e-10);
%! assert_equal (yd(13), 0.100922761449077, -1e-10);
%! assert_equal (yd(14), 0.0935406698564592, -1e-10);
%! assert_equal (yd(15), 0.0900774663932559, -1e-10);
%! assert_equal (yd(16), 0.0960013670539986, -1e-10);
%! assert_equal (yd(17), 0.118147641831852, -1e-10);
%! assert_equal (yd(18), 0.164718614718615, -1e-10);
%! assert_equal (yd(19), 0.245283663704716, -1e-10);
%! assert_equal (yd(20), 0.370779220779221, -1e-10);
%! assert_equal (get (h(2), 'YData'), [0.3, 0.3], 1e-12);
%! assert_equal (get (h(2), 'XData'), [0, n]);
%! assert_equal (get (h(2), 'LineStyle'), ':');
%! assert_equal (get (get (ax, 'xlabel'), 'string'), 'Row number');
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Leverage');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Case order plot of leverage');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'leverage', 'Color', [1 0 0]);
%! assert_equal (get (h(1), 'Color'), [1 0 0], 1e-10);
%! assert_equal (get (h(2), 'Color'), [0.8510 0.8510 0.8660], 1e-4);
%! assert_equal (get (h(2), 'LineStyle'), ':');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'cookd');
%! yd = get (h(1), 'YData');
%! assert_equal (numel (h), 2);
%! assert_equal (yd(1), 0.078517048682575, -1e-8);
%! assert_equal (yd(2), 0.077211407930332, -1e-8);
%! assert_equal (yd(3), 0.001953301452841, -1e-7);
%! assert_equal (get (h(2), 'YData'), [0.1668641787, 0.1668641787], -1e-8);
%! assert_equal (get (h(2), 'XData'), [0, n]);
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Cook''s distance');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Case order plot of Cook''s distance');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'covratio');
%! yd = get (h(1), 'YData');
%! yv = get (h(2), 'YData');
%! xv = get (h(2), 'XData');
%! assert_equal (numel (h), 2);
%! assert_equal (yd(1), 1.774933177, -1e-8);
%! assert_equal (yd(2), 1.397661919, -1e-8);
%! assert_equal (yd(3), 1.428481535, -1e-8);
%! assert_equal (numel (xv), 5);
%! assert_equal (sum (isnan (xv)), 1);
%! assert_equal (yv(1), 0.55, 1e-12);
%! assert_equal (yv(2), 0.55, 1e-12);
%! assert_equal (yv(4), 1.45, 1e-12);
%! assert_equal (yv(5), 1.45, 1e-12);
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Covariance ratio');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Case order plot of covariance ratio');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'dfbetas');
%! p = mdl.NumCoefficients;
%! yv = get (h(p+1), 'YData');
%! xv = get (h(p+1), 'XData');
%! assert_equal (numel (h), p + 1);
%! assert_equal (numel (get (h(1), 'YData')), n);
%! assert_equal (numel (get (h(2), 'YData')), n);
%! assert_equal (numel (get (h(3), 'YData')), n);
%! assert_equal (numel (xv), 5);
%! assert_equal (sum (isnan (xv)), 1);
%! assert_equal (yv(1), -0.6708203932, -1e-8);
%! assert_equal (yv(end), 0.6708203932, -1e-8);
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Scaled change in coefficients');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Case order plot of scaled change in coefficients');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'dfbetas', 'Color', [1 0 0]);
%! p = mdl.NumCoefficients;
%! for k = 1:p
%! assert_equal (get (h(k), 'Color'), [1 0 0], 1e-10);
%! endfor
%! assert_equal (get (h(p+1), 'Color'), [0.8510 0.8510 0.8660], 1e-4);
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'dffits');
%! yd = get (h(1), 'YData');
%! yv = get (h(2), 'YData');
%! xv = get (h(2), 'XData');
%! assert_equal (numel (h), 2);
%! assert_equal (yd(1), 0.476480465355394, -1e-8);
%! assert_equal (yd(2), 0.477020506700835, -1e-8);
%! assert_equal (yd(3), -0.074329411030064, -1e-7);
%! assert_equal (sum (isnan (xv)), 1);
%! assert_equal (yv(1), -0.7745966692, -1e-8);
%! assert_equal (yv(end), 0.7745966692, -1e-8);
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Scaled change in fit');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Case order plot of scaled change in fit');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 's2_i');
%! yd = get (h(1), 'YData');
%! assert_equal (numel (h), 2);
%! assert_equal (yd(1), 0.02359100986, -1e-8);
%! assert_equal (yd(2), 0.02314622330, -1e-8);
%! assert_equal (yd(3), 0.02411685408, -1e-8);
%! assert_equal (get (h(2), 'YData'), [0.02273796067, 0.02273796067], -1e-8);
%! assert_equal (get (h(2), 'XData'), [0, n]);
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Leave-one-out variance');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Case order plot of leave-one-out variance');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'contour');
%! yd = get (h(1), 'YData');
%! xd = get (h(1), 'XData');
%! assert_equal (numel (h), 2);
%! assert_equal (xd(1), 0.3707792208, -1e-8);
%! assert_equal (xd(2), 0.2452836637, -1e-8);
%! assert_equal (yd(1), 0.07562471113, -1e-8);
%! assert_equal (yd(2), 0.11059272450, -1e-8);
%! assert_equal (get (h(1), 'LineStyle'), 'none');
%! assert_equal (get (h(1), 'Marker'), 'x');
%! assert_equal (isgraphics (h(2)), true);
%! assert_equal (get (get (ax, 'xlabel'), 'string'), 'Leverage');
%! assert_equal (get (get (ax, 'ylabel'), 'string'), 'Residual');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Cook''s distance factorization');
%! close (fig);

%!test
%! me = fitlm (X, y, 'Exclude', [2, 7]);
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h_ex = plotDiagnostics (ax, me, 'cookd');
%! h_un = plotDiagnostics (ax, mdl, 'cookd');
%! ref_ex = get (h_ex(2), 'YData');
%! ref_un = get (h_un(2), 'YData');
%! assert_equal (! isequal (ref_ex, ref_un), true);
%! assert_equal (ref_ex(1), 3 * mean (me.Diagnostics.CooksDistance, 'omitnan'), 1e-12);
%! close (fig);

%!test
%! mw = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! hw = plotDiagnostics (ax, mw, 'leverage');
%! hu = plotDiagnostics (ax, mdl, 'leverage');
%! ydw = get (hw(1), 'YData');
%! ydu = get (hu(1), 'YData');
%! assert_equal (ydw(1) != ydu(1), true);
%! assert_equal (! isequal (ydw, ydu), true);
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! h = plotDiagnostics (mdl);
%! assert_equal (isequal (get (h(1), 'Parent'), gca ()), true);
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotEffects (ax, mdl);
%! xd1 = get (h(1), 'XData');
%! yd1 = get (h(1), 'YData');
%! xd2 = get (h(2), 'XData');
%! yd2 = get (h(2), 'YData');
%! xd3 = get (h(3), 'XData');
%! yd3 = get (h(3), 'YData');
%! ytl = get (ax, 'YTickLabel');
%! assert_equal (numel (h), 3);
%! assert_equal (xd1(1), 2.38302891604232, -1e-10);
%! assert_equal (xd1(2), -19.5277648300125, -1e-10);
%! assert_equal (yd1, [1 2]);
%! assert_equal (xd2(1), 1.39673712385796, -1e-10);
%! assert_equal (xd2(2), 3.36932070822668, -1e-10);
%! assert_equal (yd2, [1 1]);
%! assert_equal (xd3(1), -20.4857975891918, -1e-10);
%! assert_equal (xd3(2), -18.5697320708331, -1e-10);
%! assert_equal (yd3, [2 2]);
%! assert_equal (get (h(1), 'Color'), [0.1490 0.5490 0.8660], 1e-4);
%! assert_equal (get (h(2), 'Color'), [0.1490 0.5490 0.8660], 1e-4);
%! assert_equal (get (h(3), 'Color'), [0.1490 0.5490 0.8660], 1e-4);
%! assert_equal (get (h(1), 'Marker'), 'o');
%! assert_equal (get (h(1), 'LineStyle'), 'none');
%! assert_equal (get (h(2), 'LineStyle'), '-');
%! assert_equal (get (h(2), 'Marker'), 'none');
%! assert_equal (get (h(3), 'LineStyle'), '-');
%! assert_equal (get (h(3), 'Marker'), 'none');
%! assert_equal (mean (xd2), xd1(1), 1e-10);
%! assert_equal (mean (xd3), xd1(2), 1e-10);
%! assert_equal (get (get (ax, 'xlabel'), 'string'), 'Main Effect');
%! assert_equal (get (get (ax, 'ylabel'), 'string'), '');
%! assert_equal (get (get (ax, 'title'), 'string'), 'Main Effects Plot');
%! assert_equal (get (ax, 'YTick'), [1 2]);
%! assert_equal (ytl{1}, 'x1: 0.05 to 1');
%! assert_equal (ytl{2}, 'x2: 0.05 to 20');
%! close (fig);

%!test
%! ## 3-predictor model
%! X3 = [X, sin((1:n)' * pi / n)];
%! y3 = X3 * [3; -1; 2] + 0.1 * cos ((1:n)' * pi / 7);
%! m3 = fitlm (X3, y3);
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotEffects (ax, m3);
%! xd1 = get (h(1), 'XData');
%! yd1 = get (h(1), 'YData');
%! xd2 = get (h(2), 'XData');
%! yd2 = get (h(2), 'YData');
%! xd3 = get (h(3), 'XData');
%! yd3 = get (h(3), 'YData');
%! xd4 = get (h(4), 'XData');
%! yd4 = get (h(4), 'YData');
%! ytl = get (ax, 'YTickLabel');
%! assert_equal (numel (h), 4);
%! assert_equal (xd1(1), 8.10687671732127, -1e-10);
%! assert_equal (xd1(2), -25.4487243632125, -1e-10);
%! assert_equal (xd1(3), 0.661302203942261, -1e-10);
%! assert_equal (yd1, [1 2 3]);
%! assert_equal (xd2(1), 0.565266595687836, -1e-10);
%! assert_equal (xd2(2), 15.6484868389547, -1e-10);
%! assert_equal (yd2, [1 1]);
%! assert_equal (xd3(1), -33.3368582824351, -1e-10);
%! assert_equal (xd3(2), -17.5605904439899, -1e-10);
%! assert_equal (yd3, [2 2]);
%! assert_equal (xd4(1), -1.25582490831999, -1e-10);
%! assert_equal (xd4(2), 2.57842931620451, -1e-10);
%! assert_equal (yd4, [3 3]);
%! assert_equal (get (ax, 'YTick'), [1 2 3]);
%! assert_equal (ytl{1}, 'x1: 0.05 to 1');
%! assert_equal (ytl{2}, 'x2: 0.05 to 20');
%! assert_equal (ytl{3}, 'x3: 1.22465e-16 to 1');
%! assert_equal (mean (xd2), xd1(1), 1e-10);
%! assert_equal (mean (xd3), xd1(2), 1e-10);
%! assert_equal (mean (xd4), xd1(3), 1e-10);
%! close (fig);

%!test
%! me = fitlm (X, y, 'Exclude', [2, 7]);
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotEffects (ax, me);
%! xd1 = get (h(1), 'XData');
%! yd1 = get (h(1), 'YData');
%! xd2 = get (h(2), 'XData');
%! yd2 = get (h(2), 'YData');
%! xd3 = get (h(3), 'XData');
%! yd3 = get (h(3), 'YData');
%! ytl = get (ax, 'YTickLabel');
%! assert_equal (numel (h), 3);
%! assert_equal (xd1(1), 2.50035744908398, -1e-10);
%! assert_equal (xd1(2), -19.5912988214488, -1e-10);
%! assert_equal (yd1, [1 2]);
%! assert_equal (xd2(1), 1.40421088339552, -1e-10);
%! assert_equal (xd2(2), 3.59650401477245, -1e-10);
%! assert_equal (yd2, [1 1]);
%! assert_equal (xd3(1), -20.6333076647782, -1e-10);
%! assert_equal (xd3(2), -18.5492899781194, -1e-10);
%! assert_equal (yd3, [2 2]);
%! assert_equal (ytl{1}, 'x1: 0.05 to 1');
%! assert_equal (ytl{2}, 'x2: 0.05 to 20');
%! assert_equal (mean (xd2), xd1(1), 1e-10);
%! assert_equal (mean (xd3), xd1(2), 1e-10);
%! close (fig);

%!test
%! mw = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotEffects (ax, mw);
%! xd1 = get (h(1), 'XData');
%! yd1 = get (h(1), 'YData');
%! xd2 = get (h(2), 'XData');
%! yd2 = get (h(2), 'YData');
%! xd3 = get (h(3), 'XData');
%! yd3 = get (h(3), 'YData');
%! ytl = get (ax, 'YTickLabel');
%! assert_equal (numel (h), 3);
%! assert_equal (xd1(1), 2.51587141860715, -1e-10);
%! assert_equal (xd1(2), -19.6411669663483, -1e-10);
%! assert_equal (yd1, [1 2]);
%! assert_equal (xd2(1), 1.08491557053384, -1e-10);
%! assert_equal (xd2(2), 3.94682726668046, -1e-10);
%! assert_equal (yd2, [1 1]);
%! assert_equal (xd3(1), -20.8383905241664, -1e-10);
%! assert_equal (xd3(2), -18.4439434085302, -1e-10);
%! assert_equal (yd3, [2 2]);
%! assert_equal (ytl{1}, 'x1: 0.05 to 1');
%! assert_equal (ytl{2}, 'x2: 0.05 to 20');
%! assert_equal (mean (xd2), xd1(1), 1e-10);
%! assert_equal (mean (xd3), xd1(2), 1e-10);
%! close (fig);

%!test
%! mni = fitlm (X, y, 'Intercept', false);
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotEffects (ax, mni);
%! xd1 = get (h(1), 'XData');
%! yd1 = get (h(1), 'YData');
%! xd2 = get (h(2), 'XData');
%! yd2 = get (h(2), 'YData');
%! xd3 = get (h(3), 'XData');
%! yd3 = get (h(3), 'YData');
%! ytl = get (ax, 'YTickLabel');
%! assert_equal (numel (h), 3);
%! assert_equal (xd1(1), 2.81335053251731, -1e-10);
%! assert_equal (xd1(2), -19.8951125513936, -1e-10);
%! assert_equal (yd1, [1 2]);
%! assert_equal (xd2(1), 2.36234515544211, -1e-10);
%! assert_equal (xd2(2), 3.26435590959250, -1e-10);
%! assert_equal (yd2, [1 1]);
%! assert_equal (xd3(1), -20.4919734818287, -1e-10);
%! assert_equal (xd3(2), -19.2982516209584, -1e-10);
%! assert_equal (yd3, [2 2]);
%! assert_equal (ytl{1}, 'x1: 0.05 to 1');
%! assert_equal (ytl{2}, 'x2: 0.05 to 20');
%! assert_equal (mean (xd2), xd1(1), 1e-10);
%! assert_equal (mean (xd3), xd1(2), 1e-10);
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotEffects (ax, mdl);
%! assert_equal (isequal (get (h(1), 'Parent'), ax), true);
%! assert_equal (get (h(1), 'XData'), [2.38302891604232, -19.5277648300125], -1e-10);
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! h = plotEffects (mdl);
%! assert_equal (isequal (get (h(1), 'Parent'), gca ()), true);
%! assert_equal (get (h(1), 'XData'), [2.38302891604232, -19.5277648300125], -1e-10);
%! close (fig);

%!test
%! ## numeric predictor adjusted data
%! fig = figure ('visible', 'off');
%! h = plotAdjustedResponse (mdl, 'x1');
%! assert_equal (get (h(1), 'XData'), X(:,1)', 1e-10);
%! assert_equal (get (h(1), 'YData'), ...
%!   [-6.70590752804287, -6.54551694016554, -6.5544435914624,  -6.59143572669715, ...
%!    -6.49138418414686, -6.21712299745016, -5.89359961368025, -5.69299878673044, ...
%!    -5.67643670865536, -5.73777106454765, -5.70118778736348, -5.48284370035446, ...
%!    -5.16795154710756, -4.93243578806991, -4.88118846293094, -4.95163193306634, ...
%!    -4.97125247389768, -4.81620859768203, -4.52519034621851, -4.26384784484646], 1e-8);
%! xf = get (h(2), 'XData');
%! yf = get (h(2), 'YData');
%! assert_equal (numel (xf), 100);
%! assert_equal (xf(1:5), ...
%!   [0.05, 0.0595959595959596, 0.0691919191919192, 0.0787878787878788, 0.0883838383838384], 1e-10);
%! assert_equal (yf(1:5), ...
%!   [-6.78153223917696, -6.75746124002502, -6.73339024087308, -6.70931924172113, -6.68524824256919], 1e-8);
%! assert_equal (xf(end-4:end), ...
%!   [0.961616161616162, 0.971212121212121, 0.980808080808081, 0.99040404040404, 1], 1e-10);
%! assert_equal (yf(end-4:end), ...
%!   [-4.49478731974241, -4.47071632059047, -4.44664532143853, -4.42257432228658, -4.39850332313464], 1e-8);
%! close (fig);

%!test
%! ## title and axis labels follow the standard convention
%! fig = figure ('visible', 'off');
%! plotAdjustedResponse (mdl, 'x1');
%! assert_equal (get (get (gca, 'Title'),  'String'), 'Adjusted Response Plot');
%! assert_equal (get (get (gca, 'XLabel'), 'String'), 'x1');
%! assert_equal (get (get (gca, 'YLabel'), 'String'), 'Adjusted y');
%! close (fig);

%!test
%! ## ax routing and a second predictor
%! fig = figure ('visible', 'off');
%! ax  = axes (fig);
%! h   = plotAdjustedResponse (ax, mdl, 'x2');
%! assert_equal (isequal (get (h(1), 'Parent'), ax), true);
%! assert_equal (get (h(1), 'XData'), X(:,2)', 1e-8);
%! assert_equal (get (h(1), 'YData'), ...
%!   [1.45980865498274,  1.34795136885775,  0.968893310576047, 0.463886235373945, ...
%!    -0.00196069502564064, -0.391481514261341, -0.829623669406342, -1.48857191435397, ...
%!    -2.42944244115883, -3.5460929349136, -4.66270932857441, -5.6954484453929, ...
%!    -6.72952302895603, -7.94085753971093, -9.43434401734702, -11.14740482324, ...
%!    -12.9075262328114, -14.5908667583184, -16.23611644156, -18.0089254078756], 1e-7);
%! close (fig);

%!test
%! ## name-value arguments style the data points only
%! fig = figure ('visible', 'off');
%! h = plotAdjustedResponse (mdl, 'x1', 'Marker', 's', 'MarkerSize', 10, 'Color', 'r');
%! assert_equal (get (h(1), 'Marker'), 's');
%! assert_equal (get (h(1), 'MarkerSize'), 10);
%! assert_equal (get (h(1), 'Color'), [1 0 0]);
%! assert_equal (get (h(2), 'Marker'), 'none');
%! close (fig);

%!test
%! yn = y;
%! yn(3) = NaN;
%! mn = fitlm (X, yn);
%! fig = figure ('visible', 'off');
%! h = plotAdjustedResponse (mn, 'x1');
%! xd = get (h(1), 'XData');
%! yd = get (h(1), 'YData');
%! assert_equal (isnan (xd(3)), true);
%! assert_equal (isnan (yd(3)), true);
%! assert_equal (yd([1 2 4]), [-7.04679028889336, -6.88651148335619, -6.93287739924846], 1e-8);
%! close (fig);

%!test
%! w  = mod ((1:n)', 3) + 1;
%! mw = fitlm (X, y, 'Weights', w);
%! fig = figure ('visible', 'off');
%! h = plotAdjustedResponse (mw, 'x1');
%! assert_equal (get (h(1), 'YData'), ...
%!   [-6.66184168801736, -6.5023788020353,  -6.51285162315763, -6.55200839614801, ...
%!    -6.45473995928354, -6.18388034620285, -5.86437700397912, -5.66841468650568, ...
%!    -5.65710958583715, -5.72431938706618, -5.69423002314892, -5.482998317337,   ...
%!    -5.17583701321739, -4.9486705712372,  -4.90639103108588, -4.98642075413911, ...
%!    -5.01624601581846, -4.87202532838101, -4.59244873362587, -4.34316635689238], 1e-7);
%! close (fig);

%!test
%! ## robust regression 
%! mr = fitlm (X, y, 'RobustOpts', 'on');
%! fig = figure ('visible', 'off');
%! h = plotAdjustedResponse (mr, 'x1');
%! assert_equal (get (h(1), 'YData'), ...
%!   [-6.69986109212516, -6.53959779763556, -6.54873660457867, -6.58602575771814, ...
%!    -6.48635609533107, -6.2125616510561,  -5.88958987196639, -5.68962551195529, ...
%!    -5.67378476307741, -5.7359253104254,  -5.70023308695542, -5.48286491591908, ...
%!    -5.16903354090336, -4.93466342235539, -4.88464659996458, -4.95640543510664, ...
%!    -4.97742620320314, -4.82386741651113, -4.53441911682976, -4.27473142949835], 1e-7);
%! close (fig);

%!test
%! ## numeric predictor averaged over a categorical predictor
%! wt = [3504;3693;3436;3433;3449;3672;3705;3288;3092;2500;2700;3100];
%! yr = categorical ([70;70;70;70;70;76;76;76;82;82;82;82]);
%! mg = [18;15;18;16;17;20;22;24;30;32;28;26];
%! tc = table (mg, wt, yr, 'VariableNames', {'MPG','Weight','Year'});
%! mc = fitlm (tc, 'MPG ~ Year + Weight');
%! fig = figure ('visible', 'off');
%! h = plotAdjustedResponse (mc, 'Weight');
%! assert_equal (get (h(1), 'XData'), wt', 1e-10);
%! assert_equal (get (h(1), 'YData'), ...
%!   [22.1247351073949, 19.1247351073949, 22.1247351073949, 20.1247351073949, ...
%!    21.1247351073949, 18.6102199722546, 20.6102199722546, 22.6102199722546, ...
%!    25.8864161365654, 27.8864161365654, 23.8864161365654, 21.8864161365654], 1e-8);
%! xf = get (h(2), 'XData');
%! yf = get (h(2), 'YData');
%! assert_equal (xf(1:5), ...
%!   [2500, 2512.17171717172, 2524.34343434343, 2536.51515151515, 2548.68686868687], 1e-8);
%! assert_equal (yf(1:5), ...
%!   [26.9912481948118, 26.9176291703666, 26.8440101459214, 26.7703911214761, 26.6967720970309], 1e-8);
%! close (fig);

%!test
%! ## categorical predictor evaluated per level
%! wt = [3504;3693;3436;3433;3449;3672;3705;3288;3092;2500;2700;3100];
%! yr = categorical ([70;70;70;70;70;76;76;76;82;82;82;82]);
%! mg = [18;15;18;16;17;20;22;24;30;32;28;26];
%! tc = table (mg, wt, yr, 'VariableNames', {'MPG','Weight','Year'});
%! mc = fitlm (tc, 'MPG ~ Year + Weight');
%! fig = figure ('visible', 'off');
%! h = plotAdjustedResponse (mc, 'Year');
%! assert_equal (get (h(1), 'XData'), [1 1 1 1 1 2 2 2 3 3 3 3]);
%! assert_equal (get (h(1), 'YData'), ...
%!   [19.2479799272553, 17.3911214761304, 18.8366909043795, 16.8185458004291, ...
%!    17.9153196881646, 22.2641057484776, 24.463701891932,  23.9415324428265, ...
%!    28.7560523180671, 27.1754184718549, 24.3850920685482, 24.8044392619348], 1e-7);
%! assert_equal (get (h(2), 'XData'), [1 2 3]);
%! assert_equal (get (h(2), 'YData'), [18.0419315592718, 23.5564466944121, 26.2802505301012], 1e-8);
%! assert_equal (get (gca, 'XTickLabel'), {'70'; '76'; '82'});
%! close (fig);

%!test
%! ## added variable plot for the whole model
%! fig = figure ('visible', 'off');
%! h = plotAdded (mdl);
%! assert_equal (get (get (gca, 'Title'), 'String'), 'Added Variable Plot for Whole Model');
%! assert_equal (get (get (gca, 'XLabel'), 'String'), 'Adjusted Whole Model');
%! assert_equal (get (h(1), 'XData'), ...
%!   [0.0284033824838481, 0.0204548552857604, -0.0238455815942649, -0.104497928156226, ...
%!    -0.221502184400124, -0.374858350325959, -0.564566425933731, -0.79062641122344, ...
%!    -1.05303830619508, -1.35180211084867, -1.68691782518418, -2.05838544920164, ...
%!    -2.46620498290103, -2.91037642628236, -3.39089977934563, -3.90777504209083, ...
%!    -4.46100221451797, -5.05058129662704, -5.67651228841805, -6.338795189891], 1e-7);
%! assert_equal (get (h(1), 'YData'), ...
%!   [0.268294196961578, 0.281859485365135, 0.0282240016119717, -0.351360499061586, ...
%!    -0.691784854932629, -0.955883099639786, -1.26860268025624, -1.80212835067532, ...
%!    -2.61757630295165, -3.60880422217788, -4.59999804131014, -5.50731458360009, ...
%!    -6.41596659263467, -7.50187852886103, -8.86994243196858, -10.457580663333, ...
%!    -12.0922794983759, -13.6501974493543, -15.1700245580674, -16.8174109498545], 1e-7);
%! assert_equal (get (h(2), 'DisplayName'), 'Fit: y = 2.69267*x');
%! close (fig);

%!test
%! ## added variable plot for just the intercept term
%! fig = figure ('visible', 'off');
%! h = plotAdded (mdl, 1);
%! assert_equal (get (get (gca, 'Title'), 'String'), 'Added Variable Plot for (Intercept)');
%! assert_equal (get (h(1), 'YData'), ...
%!   [-5.41993222738086, -5.40485070721962, -5.55724508427077, -5.73586360329798, ...
%!    -5.77559710257835, -5.63927961575051, -5.45185858988763, -5.38551877888305, ...
%!    -5.50137637479139, -5.6932890627053,  -5.78544277558092, -5.6939943366699,  ...
%!    -5.50415648955918, -5.3918536946959,  -5.4619779917695,  -5.65195174215564, ...
%!    -5.78926122127593, -5.75006494138741, -5.57305294428921, -5.42387535532067], 1e-7);
%! assert_equal (get (h(2), 'DisplayName'), 'Fit: y = 0.116189*x');
%! close (fig);

%!test
%! ## added variable plot for one predictor picked by index
%! fig = figure ('visible', 'off');
%! h = plotAdded (mdl, 2);
%! assert_equal (get (get (gca, 'Title'), 'String'), 'Added Variable Plot for x1');
%! assert_equal (get (h(1), 'YData'), ...
%!   [-5.90289714234183, -5.75941203626873, -5.79651449057265, -5.87295275001727, ...
%!    -5.82361765287968, -5.6113432327985,  -5.36107693684693, -5.24500351891828, ...
%!    -5.32423917106718, -5.49264157838629, -5.57439667383173, -5.48566128065516, ...
%!    -5.31164814244354, -5.22828171964398, -5.34045405194592, -5.58558750072505, ...
%!    -5.79116834140295, -5.83335508623668, -5.75083777702536, -5.70926653910833], 1e-7);
%! assert_equal (get (h(2), 'DisplayName'), 'Fit: y = 2.50845*x');
%! xf = get (h(2), 'XData');
%! yf = get (h(2), 'YData');
%! assert_equal (xf(1:3), [0.370121951219512, 0.372449462532903, 0.374776973846294], 1e-9);
%! assert_equal (yf(1:3), [-5.97852185347592, -5.97268340425253, -5.96684495502913], 1e-7);
%! close (fig);

%!test
%! ## added variable plot for the other predictor, a negative slope this time
%! fig = figure ('visible', 'off');
%! h = plotAdded (mdl, 3);
%! assert_equal (get (get (gca, 'Title'), 'String'), 'Added Variable Plot for x2');
%! assert_equal (get (h(1), 'YData'), ...
%!   [-8.3040737600235,  -7.38815394983204, -6.7394349117973,  -6.21666489068295, ...
%!    -5.65473472476609, -5.01647844768535, -4.4268435065139,  -4.05801465514508, ...
%!    -3.9711080856335,  -4.05998148307182, -4.14882078041619, -4.15378280091823, ...
%!    -4.16008028816491, -4.34363770260337, -4.80934708392301, -5.49463079349955, ...
%!    -6.22697510675454, -6.88253853594507, -7.50001112287024, -8.2450429928694], 1e-7);
%! assert_equal (get (h(2), 'DisplayName'), 'Fit: y = -0.978835*x');
%! close (fig);

%!test
%! ## ax argument sends the plot to that axes
%! fig = figure ('visible', 'off');
%! ax  = axes (fig);
%! h   = plotAdded (ax, mdl, 2);
%! assert_equal (isequal (get (h(1), 'Parent'), ax), true);
%! close (fig);

%!test
%! ## name-value styling only changes the data points, not the fit line
%! fig = figure ('visible', 'off');
%! h = plotAdded (mdl, 2, 'Marker', 's', 'MarkerSize', 10, 'Color', 'r');
%! assert_equal (get (h(1), 'Marker'), 's');
%! assert_equal (get (h(1), 'MarkerSize'), 10);
%! assert_equal (get (h(1), 'Color'), [1 0 0]);
%! assert_equal (get (h(2), 'Marker'), 'none');
%! close (fig);

%!test
%! ## a missing observation leaves a gap in the adjusted data
%! yn = y;
%! yn(3) = NaN;
%! mn = fitlm (X, yn);
%! fig = figure ('visible', 'off');
%! h = plotAdded (mn, 2);
%! xd = get (h(1), 'XData');
%! yd = get (h(1), 'YData');
%! assert_equal (isnan (xd(3)), true);
%! assert_equal (isnan (yd(3)), true);
%! close (fig);

%!test
%! ## weighted fit still gives a proper added variable plot
%! w  = mod ((1:n)', 3) + 1;
%! mw = fitlm (X, y, 'Weights', w);
%! fig = figure ('visible', 'off');
%! h = plotAdded (mw, 2);
%! assert_equal (get (h(1), 'YData'), ...
%!   [-6.04899282212585, -5.90562605001686, -5.9429257275943,  -6.01964009962184, ...
%!    -5.97066000437658, -5.75881947549715, -5.50906596005672, -5.39358421194863, ...
%!    -5.4734904232275,  -5.64264227898597, -5.7252257121802,  -5.63739754606181, ...
%!    -5.46437052421778, -5.38206910709522, -5.49538533438357, -5.74174156745852, ...
%!    -5.94862408174164, -5.99219138949,    -5.91113353250271, -5.87110063611913], 1e-7);
%! assert_equal (get (h(2), 'DisplayName'), 'Fit: y = 2.38836*x');
%! close (fig);

%!test
%! ## robust fit still gives a proper added variable plot
%! mr = fitlm (X, y, 'RobustOpts', 'on');
%! fig = figure ('visible', 'off');
%! h = plotAdded (mr, 2);
%! assert_equal (get (h(1), 'YData'), ...
%!   [-5.89409568732029, -5.75060102429134, -5.78768755033552, -5.86410351021651, ...
%!    -5.81473974221138, -5.60243027995878, -5.35212257053188, -5.23600136782402, ...
%!    -5.31518286388981, -5.4835247438219,  -5.56521294057644, -5.47640427740507, ...
%!    -5.30231149789475, -5.2188590624926,  -5.33093901088806, -5.5759737044568,  ...
%!    -5.78144941862042, -5.82352466563597, -5.74088948730258, -5.69919400895959], 1e-7);
%! assert_equal (get (h(2), 'DisplayName'), 'Fit: y = 2.48815*x');
%! close (fig);

%!test
%! ## added variable plot for the whole model
%! load carsmall
%! Year = categorical (Model_Year);
%! tbl  = table (MPG, Weight, Year);
%! mdl1  = fitlm (tbl, 'MPG ~ Year + Weight^2');
%! fig  = figure ('visible', 'off');
%! h    = plotAdded (mdl1);
%! assert_equal (get (get (gca, 'Title'),  'String'), 'Added Variable Plot for Whole Model');
%! assert_equal (get (get (gca, 'XLabel'), 'String'), 'Adjusted Whole Model');
%! xd = get (h(1), 'XData');
%! yd = get (h(1), 'YData');
%! assert_equal (xd(1:5), ...
%!   [-4.54006581304461, -4.65629283432076, -4.49502736854252, ...
%!    -4.49300111649177, -4.50376945388336], 1e-7);
%! assert_equal (yd(1:5), [18, 15, 18, 16, 17], 1e-8);
%! assert_equal (get (h(2), 'DisplayName'), 'Fit: y = 8.44866*x');
%! xf = get (h(2), 'XData');
%! yf = get (h(2), 'YData');
%! assert_equal (xf(1:3), [-5.06005142297709, -5.03050027197254, -5.000949120968], 1e-7);
%! assert_equal (yf(1:3), [11.4556384009199, 11.7053059986795, 11.9549735964392], 1e-7);
%! close (fig);

%!test
%! ## added variable plot for the weight terms picked as a pair
%! load carsmall
%! Year = categorical (Model_Year);
%! tbl  = table (MPG, Weight, Year);
%! mdl1  = fitlm (tbl, 'MPG ~ Year + Weight^2');
%! fig  = figure ('visible', 'off');
%! h    = plotAdded (mdl1, [2 5]);
%! assert_equal (get (get (gca, 'Title'),  'String'), 'Added Variable Plot for Specified Terms');
%! assert_equal (get (get (gca, 'XLabel'), 'String'), 'Adjusted Specified Terms');
%! xd = get (h(1), 'XData');
%! yd = get (h(1), 'YData');
%! assert_equal (xd(1:5), ...
%!   [-2181.48546848357, -2241.34798193195, -2158.28850051435, ...
%!    -2157.24488312395, -2162.79109548355], 1e-6);
%! assert_equal (yd(1:5), ...
%!   [24.0284299339692, 21.0284299339692, 24.0284299339692, ...
%!    22.0284299339692, 23.0284299339692], 1e-7);
%! assert_equal (get (h(2), 'DisplayName'), 'Fit: y = 0.0164036*x');
%! close (fig);

%!test
%! ## multiple predictors delegate to plotAdded
%! fig = figure ('visible', 'off');
%! ax  = axes (fig);
%! h1  = plot (ax, mdl);
%! h2  = plotAdded (ax, mdl);
%! assert_equal (get (h1(1), 'XData'), get (h2(1), 'XData'));
%! assert_equal (get (h1(1), 'YData'), get (h2(1), 'YData'));
%! assert_equal (get (h1(2), 'YData'), get (h2(2), 'YData'));
%! assert_equal (get (h1(3), 'YData'), get (h2(3), 'YData'));
%! assert_equal (get (get (ax, 'Title'), 'String'), 'Added Variable Plot for Whole Model');
%! close (fig);

%!test
%! ## no axes argument uses the current axes
%! fig = figure ('visible', 'off');
%! h = plot (mdl);
%! assert_equal (isequal (get (h(1), 'Parent'), gca ()), true);
%! close (fig);

%!test
%! ## no predictors delegate to plotResiduals
%! mc  = fitlm (X, y, 'constant');
%! fig = figure ('visible', 'off');
%! ax  = axes (fig);
%! h1  = plot (ax, mc);
%! h2  = plotResiduals (ax, mc);
%! assert_equal (get (h1(1), 'type'), 'patch');
%! assert_equal (get (h1(1), 'YData'), get (h2(1), 'YData'));
%! assert_equal (get (get (ax, 'Title'), 'String'), 'Histogram of residuals');
%! close (fig);

%!test
%! ## single predictor: data, fit line, and confidence bounds
%! m1  = fitlm (X(:,1), y);
%! fig = figure ('visible', 'off');
%! ax  = axes (fig);
%! h   = plot (ax, m1);
%! assert_equal (numel (h), 3);
%! assert_equal (get (h(1), 'XData'), X(:,1)', 1e-15);
%! assert_equal (get (h(1), 'YData'), y', 1e-15);
%! xf = get (h(2), 'XData');
%! yf = get (h(2), 'YData');
%! assert_equal (numel (xf), 100);
%! assert_equal (xf(1:3),       [0.05, 0.0595959595959596, 0.0691919191919192], 1e-12);
%! assert_equal (xf(end-2:end), [0.980808080808081, 0.99040404040404, 1],         1e-12);
%! assert_equal (yf(1:3),       [2.98235017582927, 2.80917102518311, 2.63599187453694],   1e-7);
%! assert_equal (yf(end-2:end), [-13.8160274368485, -13.9892065874947, -14.1623857381409], 1e-7);
%! yb = get (h(3), 'YData');
%! assert_equal (numel (yb), 201);
%! assert_equal (isnan (yb(101)), true);
%! assert_equal (yb(1:3),       [1.59215527128182, 1.43944293975561, 1.28661387851973], 1e-7);
%! assert_equal (yb(102:104),   [4.37254508037672, 4.1788991106106,  3.98536987055415], 1e-7);
%! assert_equal (yb(end-2:end), [-12.4666494408313, -12.6194785020672, -12.7721908335934], 1e-7);
%! assert_equal (yb(1:100) + yb(102:201), 2 * yf, 1e-10);
%! [yp, ~] = predict (m1, xf');
%! assert_equal (yf', yp, 1e-10);
%! assert_equal (get (get (ax, 'Title'),  'String'), 'y vs. x1');
%! assert_equal (get (get (ax, 'XLabel'), 'String'), 'x1');
%! assert_equal (get (get (ax, 'YLabel'), 'String'), 'y');
%! close (fig);

%!test
%! ## real dataset with missing rows leaves gaps in the data
%! load carsmall
%! tbl2 = table (MPG, Weight);
%! mdl2 = fitlm (tbl2, 'MPG ~ Weight');
%! fig  = figure ('visible', 'off');
%! ax   = axes (fig);
%! h    = plot (ax, mdl2);
%! xd = get (h(1), 'XData');
%! yd = get (h(1), 'YData');
%! assert_equal (numel (xd), 100);
%! assert_equal (any (isnan (xd)), true);
%! assert_equal (any (isnan (yd)), true);
%! assert_equal (xd(1:3),       [3504, 3693, 3436], 1e-10);
%! assert_equal (xd(end-2:end), [2295, 2625, 2720], 1e-10);
%! assert_equal (yd(1:3),       [18, 15, 18], 1e-10);
%! assert_equal (yd(end-2:end), [32, 28, 31], 1e-10);
%! xf = get (h(2), 'XData');
%! yf = get (h(2), 'YData');
%! assert_equal (numel (xf), 100);
%! assert_equal (xf(1:3),       [1795, 1824.666667, 1854.333333], 1e-4);
%! assert_equal (xf(end-2:end), [4672.666667, 4702.333333, 4732], 1e-4);
%! assert_equal (yf(1:3),       [33.77920696, 33.52371956, 33.26823216],   1e-6);
%! assert_equal (yf(end-2:end), [8.996929296, 8.741441898, 8.485954499],   1e-6);
%! yb = get (h(3), 'YData');
%! assert_equal (numel (yb), 201);
%! assert_equal (yb(1:3),       [32.2768265, 32.0472587, 31.81746955],    1e-6);
%! assert_equal (yb(end-2:end), [11.0004001, 10.77351303, 10.54671087],   1e-6);
%! assert_equal (get (get (ax, 'Title'),  'String'), 'MPG vs. Weight');
%! assert_equal (get (get (ax, 'XLabel'), 'String'), 'Weight');
%! assert_equal (get (get (ax, 'YLabel'), 'String'), 'MPG');
%! close (fig);

%!test
%! ## categorical predictor: group codes with per-level fit and bounds
%! yv = [4.73087805313537; 7.43361607881479; 4.48799323173627; 5.16869512618961; ...
%!       6.50195295169252; 3.73839298415678; 4.3512762614163;  7.45869429457297; ...
%!       4.08830081430877; 5.37758996783874; 6.70425002011403; 4.92219481850025; ...
%!       5.90846112458998; 6.93808332486038; 3.44469932246968; 4.65954705986638; ...
%!       7.00708466319882; 3.97022469477047; 4.66949449276983; 7.15297545753449; ...
%!       3.79547107710474; 4.35907258855759; 6.85754858526846; 3.96760657205439; ...
%!       5.50019159441948; 6.59933535084439; 4.04590000762589; 5.13690239596885; ...
%!       5.93326622029102; 4.20198864027471];
%! grp  = categorical (repmat ({'A';'B';'C'}, 10, 1));
%! tbl4 = table (yv, grp, 'VariableNames', {'Response','Group'});
%! mdl4 = fitlm (tbl4, 'Response ~ Group');
%! fig  = figure ('visible', 'off');
%! ax   = axes (fig);
%! h    = plot (ax, mdl4);
%! assert_equal (numel (h), 3);
%! assert_equal (get (h(1), 'XData'), repmat ([1 2 3], 1, 10));
%! assert_equal (get (h(1), 'YData'), yv', 1e-14);
%! assert_equal (get (h(2), 'XData'), [1 2 3]);
%! assert_equal (get (h(2), 'YData'), [4.98621086647521, 6.85868069471919, 4.0662772163002], 1e-10);
%! assert_equal (get (h(3), 'XData'), [1 2 3 NaN 1 2 3]);
%! assert_equal (get (h(3), 'YData'), ...
%!   [4.68577486025387, 6.55824468849784, 3.76584121007885, NaN, ...
%!    5.28664687269656, 7.15911670094053, 4.36671322252154], 1e-10);
%! assert_equal (get (ax, 'XTick'), [1 2 3]);
%! assert_equal (get (ax, 'XTickLabel'), {'A'; 'B'; 'C'});
%! assert_equal (get (get (ax, 'Title'),  'String'), 'Response vs. Group');
%! assert_equal (get (get (ax, 'XLabel'), 'String'), 'Group');
%! assert_equal (get (get (ax, 'YLabel'), 'String'), 'Response');
%! close (fig);

%!test
%! ## continuous by continuous, effects mode
%! mi  = fitlm (X, y, 'y ~ x1*x2');
%! fig = figure ('visible', 'off');
%! h   = plotInteraction (mi, 'x1', 'x2');
%! assert_equal (numel (h), 11);
%! assert_equal (get (h(1), 'XData'), [1.76843380852813, -18.8740348676059], 1e-9);
%! assert_equal (get (h(1), 'YData'), [1, 7]);
%! assert_equal (get (h(2), 'XData'), [-2.25617028020123, 5.79303789725749], 1e-9);
%! assert_equal (get (h(3), 'XData'), [-23.1321080641968, -14.615961671015], 1e-9);
%! assert_equal (get (h(4), 'XData'), ...
%!   [1.97967308209488, 1.68393809910143, 1.38820311610799], 1e-9);
%! assert_equal (get (h(4), 'YData'), [2, 3, 4]);
%! assert_equal (get (h(5), 'XData'), [-0.771057026434185, 4.73040319062394], 1e-9);
%! assert_equal (get (h(6), 'XData'), [-2.86059582662229, 6.22847202482516], 1e-9);
%! assert_equal (get (h(7), 'XData'), [-4.99614754441031, 7.77255377662629], 1e-9);
%! assert_equal (get (h(8), 'XData'), ...
%!   [-18.5782998846125, -18.8740348676059, -19.1697698505994], 1e-9);
%! assert_equal (get (h(8), 'YData'), [8, 9, 10]);
%! assert_equal (get (h(9), 'XData'), [-24.674318784563, -12.4822809846619], 1e-9);
%! assert_equal (get (h(10), 'XData'), [-23.1321080641968, -14.615961671015], 1e-9);
%! assert_equal (get (h(11), 'XData'), [-21.6439971135684, -16.6955425876303], 1e-9);
%! assert_equal (get (h(1), 'Tag'), 'main');
%! assert_equal (get (h(4), 'Tag'), 'conditional1');
%! assert_equal (get (h(8), 'Tag'), 'conditional2');
%! ax = gca ();
%! assert_equal (get (get (ax, 'Title'), 'String'), 'Interaction of x1 and x2');
%! assert_equal (get (get (ax, 'XLabel'), 'String'), 'Effect');
%! assert_equal (get (ax, 'YTick'), [1, 2, 3, 4, 7, 8, 9, 10]);
%! assert_equal (get (ax, 'YTickLabel'), ...
%!   {'x1: 0.05 to 1'; 'x2=0.05'; 'x2=10.025'; 'x2=20'; ...
%!    'x2: 0.05 to 20'; 'x1=0.05'; 'x1=0.525'; 'x1=1'});
%! assert_equal (get (ax, 'YLim'), [0.5, 10.5]);
%! close (fig);

%!test
%! ## continuous by continuous, predictions mode
%! mi  = fitlm (X, y, 'y ~ x1*x2');
%! fig = figure ('visible', 'off');
%! h   = plotInteraction (mi, 'x1', 'x2', 'predictions');
%! assert_equal (numel (h), 3);
%! xd = get (h(1), 'XData');
%! assert_equal (numel (xd), 101);
%! assert_equal (xd(1:3), [0.05, 0.2495, 0.449], 1e-9);
%! assert_equal (xd(end-2:end), [19.601, 19.8005, 20], 1e-9);
%! yd1 = get (h(1), 'YData');
%! assert_equal (yd1(1:3), ...
%!   [0.215349913094656, 0.0295669142485309, -0.156216084597594], 1e-9);
%! assert_equal (yd1(end-2:end), ...
%!   [-17.9913839738256, -18.1771669726717, -18.3629499715178], 1e-9);
%! yd2 = get (h(2), 'YData');
%! assert_equal (yd2(1:3), [1.20518645414209, 1.01644610546604, 0.827705756789976], 1e-9);
%! assert_equal (yd2(end-2:end), ...
%!   [-17.2913677161117, -17.4801080647878, -17.6688484134638], 1e-9);
%! yd3 = get (h(3), 'YData');
%! assert_equal (yd3(1:3), [2.19502299518953, 2.00332529668354, 1.81162759817755], 1e-9);
%! assert_equal (yd3(end-2:end), ...
%!   [-16.5913514583978, -16.7830491569038, -16.9747468554098], 1e-9);
%! assert_equal (get (h(1), 'DisplayName'), '0.05');
%! assert_equal (get (h(2), 'DisplayName'), '0.525');
%! assert_equal (get (h(3), 'DisplayName'), '1');
%! ax = gca ();
%! assert_equal (get (get (ax, 'Title'), 'String'), 'Interaction of x1 and x2');
%! assert_equal (get (get (ax, 'XLabel'), 'String'), 'x2');
%! assert_equal (get (get (ax, 'YLabel'), 'String'), 'Adjusted y');
%! close (fig);

%!test
%! ## swapping var1/var2 order swaps roles and title
%! mi  = fitlm (X, y, 'y ~ x1*x2');
%! fig = figure ('visible', 'off');
%! h   = plotInteraction (mi, 'x2', 'x1');
%! assert_equal (numel (h), 11);
%! assert_equal (get (h(1), 'XData'), [-18.8740348676059, 1.76843380852813], 1e-9);
%! assert_equal (get (h(4), 'XData'), ...
%!   [-18.5782998846125, -18.8740348676059, -19.1697698505994], 1e-9);
%! assert_equal (get (h(8), 'XData'), ...
%!   [1.97967308209488, 1.68393809910143, 1.38820311610799], 1e-9);
%! ax = gca ();
%! assert_equal (get (get (ax, 'Title'), 'String'), 'Interaction of x2 and x1');
%! assert_equal (get (ax, 'YTickLabel'), ...
%!   {'x2: 0.05 to 20'; 'x1=0.05'; 'x1=0.525'; 'x1=1'; ...
%!    'x1: 0.05 to 1'; 'x2=0.05'; 'x2=10.025'; 'x2=20'});
%! close (fig);

%!test
%! ## interaction effects: variables given as indices into VariableNames
%! mi  = fitlm (X, y, 'y ~ x1*x2');
%! fig = figure ('visible', 'off');
%! h   = plotInteraction (mi, 1, 2);
%! assert_equal (numel (h), 11);
%! assert_equal (get (h(1), 'XData'), [1.76843380852813, -18.8740348676059], 1e-9);
%! assert_equal (get (h(1), 'YData'), [1, 7]);
%! ax = gca ();
%! assert_equal (get (get (ax, 'Title'), 'String'), 'Interaction of x1 and x2');
%! close (fig);

%!test
%! ## interaction effects: explicit axes argument is honored
%! mi  = fitlm (X, y, 'y ~ x1*x2');
%! fig = figure ('visible', 'off');
%! axtarget = axes (fig);
%! h   = plotInteraction (axtarget, mi, 'x1', 'x2');
%! assert_equal (numel (h), 11);
%! assert_equal (isequal (get (h(1), 'Parent'), axtarget), true);
%! assert_equal (isequal (gca (), axtarget), true);
%! close (fig);

%!test
%! ## no interaction term: conditional effects collapse to the main effect
%! mn  = fitlm (X, y, 'y ~ x1 + x2');
%! fig = figure ('visible', 'off');
%! h   = plotInteraction (mn, 'x1', 'x2');
%! xd1  = get (h(1), 'XData');
%! eff1 = xd1(1);
%! eff2 = xd1(2);
%! assert_equal (eff1, 2.38302891604232, 1e-9);
%! assert_equal (eff2, -19.5277648300125, 1e-9);
%! assert_equal (get (h(4), 'XData'), [eff1, eff1, eff1], 1e-9);
%! assert_equal (get (h(8), 'XData'), [eff2, eff2, eff2], 1e-9);
%! close (fig);

%!test
%! ## categorical by continuous, effects mode
%! xc   = (1:30)' / 30;
%! grp  = categorical (repmat ({'A';'B';'C'}, 10, 1));
%! yv   = 2*xc + 3*double (grp == 'B') - 1*double (grp == 'C') + ...
%!        1.5*xc.*double (grp == 'B') + 0.3*sin ((1:30)');
%! tblc = table (yv, xc, grp, 'VariableNames', {'Response','Xc','Group'});
%! mdlc = fitlm (tblc, 'Response ~ Xc*Group');
%! fig  = figure ('visible', 'off');
%! h    = plotInteraction (mdlc, 'Group', 'Xc');
%! assert_equal (numel (h), 11);
%! assert_equal (get (h(1), 'XData'), [4.7896970899464, 2.32528157787528], 1e-9);
%! assert_equal (get (h(2), 'XData'), [4.56862113373247, 5.01077304616034], 1e-9);
%! assert_equal (get (h(3), 'XData'), [2.02254960835685, 2.62801354739371], 1e-9);
%! assert_equal (get (h(4), 'XData'), ...
%!   [4.08328517685389, 4.7896970899464, 5.49610900303892], 1e-9);
%! assert_equal (get (h(5), 'XData'), [3.64076372661607, 4.52580662709171], 1e-9);
%! assert_equal (get (h(6), 'XData'), [4.56862113373247, 5.01077304616034], 1e-9);
%! assert_equal (get (h(7), 'XData'), [5.07555715247022, 5.91666085360761], 1e-9);
%! assert_equal (get (h(8), 'XData'), ...
%!   [1.88553612240401, 3.25156621870343, 1.8387423925184], 1e-9);
%! assert_equal (get (h(9), 'XData'), [1.36118897012269, 2.40988327468532], 1e-9);
%! assert_equal (get (h(10), 'XData'), [2.72721906642211, 3.77591337098474], 1e-9);
%! assert_equal (get (h(11), 'XData'), [1.31439524023708, 2.36308954479972], 1e-9);
%! ax = gca ();
%! assert_equal (get (get (ax, 'Title'), 'String'), 'Interaction of Group and Xc');
%! assert_equal (get (ax, 'YTick'), [1, 2, 3, 4, 7, 8, 9, 10]);
%! assert_equal (get (ax, 'YTickLabel'), ...
%!   {'Group: C to B'; 'Xc=0.0333333'; 'Xc=0.516667'; 'Xc=1'; ...
%!    'Xc: 0.033333 to 1'; 'Group=A'; 'Group=B'; 'Group=C'});
%! close (fig);

%!test
%! ## categorical by continuous, predictions mode
%! xc   = (1:30)' / 30;
%! grp  = categorical (repmat ({'A';'B';'C'}, 10, 1));
%! yv   = 2*xc + 3*double (grp == 'B') - 1*double (grp == 'C') + ...
%!        1.5*xc.*double (grp == 'B') + 0.3*sin ((1:30)');
%! tblc = table (yv, xc, grp, 'VariableNames', {'Response','Xc','Group'});
%! mdlc = fitlm (tblc, 'Response ~ Xc*Group');
%! fig  = figure ('visible', 'off');
%! h    = plotInteraction (mdlc, 'Group', 'Xc', 'predictions');
%! assert_equal (numel (h), 3);
%! xd = get (h(1), 'XData');
%! assert_equal (numel (xd), 101);
%! assert_equal (xd(1:3), [0.0333333333333333, 0.043, 0.0526666666666667], 1e-9);
%! assert_equal (xd(end-2:end), [0.980666666666667, 0.990333333333333, 1], 1e-9);
%! yd1 = get (h(1), 'YData');
%! assert_equal (yd1(1:3), [0.107201421526318, 0.126056782750359, 0.144912143974399], 1e-9);
%! assert_equal (yd1(end-2:end), [1.95502682148225, 1.97388218270629, 1.99273754393033], 1e-9);
%! yd2 = get (h(2), 'YData');
%! assert_equal (yd2(1:3), [3.18658823804771, 3.21910390023474, 3.25161956242178], 1e-9);
%! assert_equal (yd2(end-2:end), [6.37312313237707, 6.4056387945641, 6.43815445675114], 1e-9);
%! yd3 = get (h(3), 'YData');
%! assert_equal (yd3(1:3), [-0.896696938806178, -0.878309514880994, -0.85992209095581], 1e-9);
%! assert_equal (yd3(end-2:end), [0.905270605861853, 0.923658029787037, 0.942045453712221], 1e-9);
%! assert_equal (get (h(1), 'DisplayName'), 'A');
%! assert_equal (get (h(2), 'DisplayName'), 'B');
%! assert_equal (get (h(3), 'DisplayName'), 'C');
%! ax = gca ();
%! assert_equal (get (get (ax, 'Title'), 'String'), 'Interaction of Group and Xc');
%! assert_equal (get (get (ax, 'XLabel'), 'String'), 'Xc');
%! assert_equal (get (get (ax, 'YLabel'), 'String'), 'Adjusted Response');
%! close (fig);

%!test
%! ## default h on shared fixture, no polynomial hierarchy present
%! t = anova (mdl);
%! assert_equal (t.Properties.RowNames, {'x1'; 'x2'; 'Error'});
%! assert_equal (t.SumSq(1), 0.590865029026992, -1e-9);
%! assert_equal (t.DF(1), 1);
%! assert_equal (t.MeanSq(1), 0.590865029026992, -1e-9);
%! assert_equal (t.F(1), 25.9858409295, -1e-8);
%! assert_equal (t.pValue(1), 8.93779416897245e-05, -1e-8);
%! assert_equal (t.SumSq(2), 42.0518254818947, -1e-9);
%! assert_equal (t.DF(2), 1);
%! assert_equal (t.MeanSq(2), 42.0518254818947, -1e-9);
%! assert_equal (t.F(2), 1849.41059985747, -1e-7);
%! assert_equal (t.pValue(2), 8.65693830575066e-19, -1e-8);
%! assert_equal (t.SumSq(3), 0.386545331386823, -1e-9);
%! assert_equal (t.DF(3), 17);
%! assert_equal (t.MeanSq(3), 0.0227379606698131, -1e-9);
%! assert (isnan (t.F(3)));
%! assert (isnan (t.pValue(3)));

%!test
%! ## explicit Type 2 on shared fixture reaches the identical result as h
%! t = anova (mdl, 'components', 2);
%! assert_equal (t.SumSq(1), 0.590865029026992, -1e-9);
%! assert_equal (t.MeanSq(1), 0.590865029026992, -1e-9);
%! assert_equal (t.F(1), 25.9858409295, -1e-8);
%! assert_equal (t.pValue(1), 8.93779416897245e-05, -1e-8);
%! assert_equal (t.SumSq(2), 42.0518254818947, -1e-9);
%! assert_equal (t.MeanSq(2), 42.0518254818947, -1e-9);
%! assert_equal (t.F(2), 1849.41059985747, -1e-7);
%! assert_equal (t.pValue(2), 8.65693830575066e-19, -1e-8);
%! assert_equal (t.SumSq(3), 0.386545331386823, -1e-9);
%! assert_equal (t.MeanSq(3), 0.0227379606698131, -1e-9);

%!test
%! ## explicit Type 1 on shared fixture, order-dependent x1 row differs from h
%! t = anova (mdl, 'components', 1);
%! assert_equal (t.SumSq(1), 541.472049189064, -1e-9);
%! assert_equal (t.MeanSq(1), 541.472049189064, -1e-9);
%! assert_equal (t.F(1), 23813.5713686901, -1e-7);
%! assert_equal (t.pValue(1), 3.41699381539186e-28, -1e-8);
%! assert_equal (t.SumSq(2), 42.0518254818947, -1e-9);
%! assert_equal (t.F(2), 1849.41059985747, -1e-7);
%! assert_equal (t.SumSq(3), 0.386545331386823, -1e-9);
%! assert_equal (t.DF(3), 17);

%!test
%! ## summary table on shared fixture, no Linear/Nonlinear split
%! t = anova (mdl, 'summary');
%! assert_equal (t.Properties.RowNames, {'Total'; 'Model'; 'Residual'});
%! assert_equal (t.SumSq(1), 583.910420002346, -1e-9);
%! assert_equal (t.DF(1), 19);
%! assert_equal (t.MeanSq(1), 30.7321273685445, -1e-9);
%! assert (isnan (t.F(1)));
%! assert_equal (t.SumSq(2), 583.523874670959, -1e-9);
%! assert_equal (t.DF(2), 2);
%! assert_equal (t.MeanSq(2), 291.761937335479, -1e-9);
%! assert_equal (t.F(2), 12831.4909842738, -1e-7);
%! assert_equal (t.pValue(2), 9.48988083209278e-28, -1e-8);
%! assert_equal (t.SumSq(3), 0.386545331386823, -1e-9);
%! assert_equal (t.DF(3), 17);
%! assert_equal (t.MeanSq(3), 0.0227379606698131, -1e-9);

%!test
%! ## polynomial hierarchy, h and Type 2 diverge only on the Age row
%! Age = [25;31;42;29;55;38;46;33;27;50;41;36;48;30;44];
%! Sex = categorical ({'M';'F';'F';'M';'M';'F';'M';'F';'F';'M';'F';'M';'F';'M';'F'});
%! BP  = [118;122;135;120;150;128;140;124;119;145;130;126;138;121;136];
%! T = table (Age, Sex, BP);
%! m = fitlm (T, 'BP ~ Sex + Age^2');
%! t  = anova (m);
%! t2 = anova (m, 'components', 2);
%! assert_equal (t.Properties.RowNames, {'Age'; 'Sex'; 'Age^2'; 'Error'});
%! assert_equal (t.SumSq(1), 1363.92365918468, -1e-9);
%! assert_equal (t.DF(1), 1);
%! assert_equal (t.MeanSq(1), 1363.92365918468, -1e-9);
%! assert_equal (t.F(1), 707.420098987802, -1e-7);
%! assert_equal (t.pValue(1), 2.46488912775244e-11, -1e-8);
%! assert_equal (t.SumSq(2), 1.90509548061236, -1e-9);
%! assert_equal (t.F(2), 0.988107233422164, -1e-9);
%! assert_equal (t.pValue(2), 0.341568882808456, -1e-9);
%! assert_equal (t.SumSq(3), 8.58235117456131, -1e-9);
%! assert_equal (t.F(3), 4.45136916320193, -1e-8);
%! assert_equal (t.pValue(3), 0.0585944334523874, -1e-9);
%! assert_equal (t.SumSq(4), 21.2082753550521, -1e-9);
%! assert_equal (t.DF(4), 11);
%! assert_equal (t.MeanSq(4), 1.92802503227746, -1e-9);
%! assert_equal (t2.SumSq(1), 0.209727728295816, -1e-8);
%! assert_equal (t2.F(1), 0.108778529731057, -1e-9);
%! assert_equal (t2.pValue(1), 0.747734511582268, -1e-9);
%! assert_equal (t2.SumSq(2), 1.90509548061236, -1e-9);
%! assert_equal (t2.SumSq(3), 8.58235117456131, -1e-9);
%! assert_equal (t2.SumSq(4), 21.2082753550521, -1e-9);

%!test
%! ## summary table with Linear/Nonlinear split, same polynomial model
%! Age = [25;31;42;29;55;38;46;33;27;50;41;36;48;30;44];
%! Sex = categorical ({'M';'F';'F';'M';'M';'F';'M';'F';'F';'M';'F';'M';'F';'M';'F'});
%! BP  = [118;122;135;120;150;128;140;124;119;145;130;126;138;121;136];
%! T = table (Age, Sex, BP);
%! m = fitlm (T, 'BP ~ Sex + Age^2');
%! t = anova (m, 'summary');
%! assert_equal (t.Properties.RowNames, ...
%!   {'Total'; 'Model'; '. Linear'; '. Nonlinear'; 'Residual'});
%! assert_equal (t.SumSq(1), 1415.73333333333, -1e-9);
%! assert_equal (t.DF(1), 14);
%! assert_equal (t.MeanSq(1), 101.12380952381, -1e-9);
%! assert_equal (t.SumSq(2), 1394.52505797828, -1e-9);
%! assert_equal (t.DF(2), 3);
%! assert_equal (t.F(2), 241.097329241452, -1e-7);
%! assert_equal (t.pValue(2), 2.58928117898032e-10, -1e-8);
%! assert_equal (t.SumSq(3), 1385.94270680372, -1e-9);
%! assert_equal (t.DF(3), 2);
%! assert_equal (t.F(3), 359.420309280577, -1e-7);
%! assert_equal (t.pValue(3), 9.54785035758658e-11, -1e-8);
%! assert_equal (t.SumSq(4), 8.58235117456131, -1e-9);
%! assert_equal (t.DF(4), 1);
%! assert_equal (t.F(4), 4.45136916320193, -1e-8);
%! assert_equal (t.pValue(4), 0.0585944334523874, -1e-9);
%! assert_equal (t.SumSq(5), 21.2082753550521, -1e-9);
%! assert_equal (t.DF(5), 11);
%! assert_equal (t.MeanSq(5), 1.92802503227746, -1e-9);

%!test
%! ## unbalanced categorical interaction, hierarchical model, default h
%! grpA = categorical ([1;1;2;1;2;1;1;2;1;1;2;1;2;2;1;1;2;1;1;2]);
%! grpB = categorical ([1;2;1;1;2;2;1;1;2;1;1;2;2;1;2;1;1;2;2;1]);
%! xv   = (1:20)' / 10;
%! yv   = 2*(grpA=='2') + 1.5*(grpB=='2') + 0.7*xv + ...
%!        0.9*(grpA=='2').*(grpB=='2') + ...
%!        [0.1;-0.2;0.05;0.15;-0.1;0.2;-0.05;0.1;0.0;-0.15; ...
%!         0.1;0.05;-0.2;0.15;0.0;-0.1;0.05;0.2;-0.05;0.1];
%! T = table (grpA, grpB, xv, yv);
%! m  = fitlm (T, 'yv ~ grpA*grpB');
%! t  = anova (m);
%! t2 = anova (m, 'components', 1);
%! assert_equal (t.Properties.RowNames, {'grpA'; 'grpB'; 'grpA:grpB'; 'Error'});
%! assert_equal (t.SumSq(1), 26.0224323327616, -1e-9);
%! assert_equal (t.F(1), 138.504723473419, -1e-8);
%! assert_equal (t.pValue(1), 2.72592668860896e-09, -1e-8);
%! assert_equal (t.SumSq(2), 15.2365308176101, -1e-9);
%! assert_equal (t.F(2), 81.0966269640544, -1e-8);
%! assert_equal (t.pValue(2), 1.15578182042777e-07, -1e-8);
%! assert_equal (t.SumSq(3), 0.0142868014375564, -1e-8);
%! assert_equal (t.F(3), 0.0760416803903897, -1e-8);
%! assert_equal (t.pValue(3), 0.786264832724038, -1e-9);
%! assert_equal (t.SumSq(4), 3.00609904761905, -1e-9);
%! assert_equal (t.DF(4), 16);
%! assert_equal (t.MeanSq(4), 0.18788119047619, -1e-9);
%! assert_equal (t2.SumSq(1), 16.3540833333333, -1e-9);
%! assert_equal (t2.F(1), 87.0448142886652, -1e-8);
%! assert_equal (t2.pValue(1), 7.14746180184104e-08, -1e-8);
%! assert_equal (t2.SumSq(2), 15.2365308176101, -1e-9);
%! assert_equal (t2.SumSq(3), 0.0142868014375564, -1e-8);

%!test
%! ## non-hierarchical model, missing grpB main effect, h still succeeds
%! grpA = categorical ([1;1;2;1;2;1;1;2;1;1;2;1;2;2;1;1;2;1;1;2]);
%! grpB = categorical ([1;2;1;1;2;2;1;1;2;1;1;2;2;1;2;1;1;2;2;1]);
%! xv   = (1:20)' / 10;
%! yv   = 2*(grpA=='2') + 1.5*(grpB=='2') + 0.7*xv + ...
%!        0.9*(grpA=='2').*(grpB=='2') + ...
%!        [0.1;-0.2;0.05;0.15;-0.1;0.2;-0.05;0.1;0.0;-0.15; ...
%!         0.1;0.05;-0.2;0.15;0.0;-0.1;0.05;0.2;-0.05;0.1];
%! T = table (grpA, grpB, xv, yv);
%! m = fitlm (T, 'yv ~ grpA + grpA:grpB');
%! t = anova (m);
%! assert_equal (t.Properties.RowNames, {'grpA'; 'grpA:grpB'; 'Error'});
%! assert_equal (t.SumSq(1), 16.3540833333333, -1e-9);
%! assert_equal (t.F(1), 22.0110535802411, -1e-8);
%! assert_equal (t.pValue(1), 0.000209917579961884, -1e-8);
%! assert_equal (t.SumSq(2), 5.62601666666667, -1e-9);
%! assert_equal (t.F(2), 7.57208776360619, -1e-8);
%! assert_equal (t.pValue(2), 0.0136181560209473, -1e-9);
%! assert_equal (t.SumSq(3), 12.6309, -1e-9);
%! assert_equal (t.DF(3), 17);
%! assert_equal (t.MeanSq(3), 0.742994117647059, -1e-9);

%!test
%! ## no-intercept model
%! grpA = categorical ([1;1;2;1;2;1;1;2;1;1;2;1;2;2;1;1;2;1;1;2]);
%! grpB = categorical ([1;2;1;1;2;2;1;1;2;1;1;2;2;1;2;1;1;2;2;1]);
%! xv   = (1:20)' / 10;
%! yv   = 2*(grpA=='2') + 1.5*(grpB=='2') + 0.7*xv + ...
%!        0.9*(grpA=='2').*(grpB=='2') + ...
%!        [0.1;-0.2;0.05;0.15;-0.1;0.2;-0.05;0.1;0.0;-0.15; ...
%!         0.1;0.05;-0.2;0.15;0.0;-0.1;0.05;0.2;-0.05;0.1];
%! T = table (grpA, grpB, xv, yv);
%! m = fitlm (T, 'yv ~ grpA + grpB - 1');
%! t = anova (m);
%! assert_equal (t.SumSq(1), 63.3745141509434, -1e-8);
%! assert_equal (t.F(1), 178.349190204051, -1e-7);
%! assert_equal (t.pValue(1), 3.91187977275845e-12, -1e-8);
%! assert_equal (t.SumSq(2), 15.2365308176101, -1e-8);
%! assert_equal (t.F(2), 85.7575941763448, -1e-7);
%! assert_equal (t.pValue(2), 4.71596005128041e-08, -1e-8);
%! assert_equal (t.SumSq(3), 3.02038584905660, -1e-8);
%! assert_equal (t.DF(3), 17);
%! assert_equal (t.MeanSq(3), 0.177669755826859, -1e-9);

%!test
%! ## weighted fit
%! grpA = categorical ([1;1;2;1;2;1;1;2;1;1;2;1;2;2;1;1;2;1;1;2]);
%! grpB = categorical ([1;2;1;1;2;2;1;1;2;1;1;2;2;1;2;1;1;2;2;1]);
%! xv   = (1:20)' / 10;
%! yv   = 2*(grpA=='2') + 1.5*(grpB=='2') + 0.7*xv + ...
%!        0.9*(grpA=='2').*(grpB=='2') + ...
%!        [0.1;-0.2;0.05;0.15;-0.1;0.2;-0.05;0.1;0.0;-0.15; ...
%!         0.1;0.05;-0.2;0.15;0.0;-0.1;0.05;0.2;-0.05;0.1];
%! T = table (grpA, grpB, xv, yv);
%! w = [1.2;0.8;1.5;1.0;0.9;1.1;1.3;0.7;1.0;1.4; ...
%!      0.8;1.2;1.0;1.1;0.9;1.3;0.7;1.5;1.0;1.2];
%! m = fitlm (T, 'yv ~ grpA + grpB', 'Weights', w);
%! t = anova (m);
%! assert_equal (t.SumSq(1), 26.6364861423312, -1e-9);
%! assert_equal (t.F(1), 134.770391630163, -1e-8);
%! assert_equal (t.pValue(1), 1.66729622333192e-09, -1e-8);
%! assert_equal (t.SumSq(2), 17.4880599694593, -1e-9);
%! assert_equal (t.F(2), 88.4828681359069, -1e-8);
%! assert_equal (t.pValue(2), 3.76674631526712e-08, -1e-8);
%! assert_equal (t.SumSq(3), 3.35993877395756, -1e-9);
%! assert_equal (t.DF(3), 17);
%! assert_equal (t.MeanSq(3), 0.197643457291621, -1e-9);

%!test
%! grp3 = categorical ([1;1;1;1;1;2;2;2;2;2;2;2;2;2;3;3;3;3;3;3]);
%! yy   = [2.1;1.9;2.3;2.0;1.8;4.1;4.3;3.9;4.0;4.2; ...
%!         4.4;3.8;4.1;4.0;6.1;5.9;6.3;6.0;5.8;6.2];
%! T = table (grp3, yy);
%! m = fitlm (T, 'yy ~ grp3');
%! t = anova (m);
%! assert_equal (t.Properties.RowNames, {'grp3'; 'Error'});
%! assert_equal (t.DF(1), 2);
%! assert_equal (t.SumSq(1), 44.3761111111111, -1e-9);
%! assert_equal (t.MeanSq(1), 22.1880555555556, -1e-9);
%! assert_equal (t.F(1), 616.446794988197, -1e-7);
%! assert_equal (t.pValue(1), 1.36582477075565e-16, -1e-8);
%! assert_equal (t.SumSq(2), 0.611888888888889, -1e-9);
%! assert_equal (t.DF(2), 17);
%! assert_equal (t.MeanSq(2), 0.0359934640522876, -1e-9);

%!test
%! ## type 3 gives different numbers than type 2 here
%! G3 = categorical ([1;1;1;1;2;2;2;3;3;3;3;1;2;3]);
%! G2 = categorical ([1;1;2;2;1;2;2;1;1;2;2;2;1;1]);
%! y1  = [10;12;15;14;9;11;13;16;18;20;22;17;10;19];
%! T1 = table (G3, G2, y1);
%! mdl1 = fitlm (T1, 'y1 ~ G3*G2');
%! t1 = anova (mdl1, 'components', 3);
%! assert_equal (t1.Properties.RowNames, {'G3'; 'G2'; 'G3:G2'; 'Error'});
%! assert_equal (t1.SumSq(1), 176.678431372549, -1e-9);
%! assert_equal (t1.F(1), 44.6345510835913, -1e-8);
%! assert_equal (t1.pValue(1), 4.57572925304662e-05, -1e-8);
%! assert_equal (t1.SumSq(2), 38.7604166666666, -1e-9);
%! assert_equal (t1.F(2), 19.5842105263158, -1e-8);
%! assert_equal (t1.pValue(2), 0.00221050293676794, -1e-9);
%! assert_equal (t1.SumSq(3), 1.85490196078431, -1e-9);
%! assert_equal (t1.SumSq(4), 15.8333333333333, -1e-9);
%! assert_equal (t1.DF(4), 8);
%! assert_equal (t1.MeanSq(4), 1.97916666666667, -1e-9);

%!test
%! ## no intercept; every level is coded, so the error row equals mdl.SSE
%! G = categorical ([1;1;1;2;2;2;3;3;3;3]);
%! y1 = [10;12;11;20;22;19;15;17;16;14];
%! T = table (G, y1);
%! mdl1 = fitlm (T, 'y1 ~ G - 1');
%! t = anova (mdl1, 'components', 3);
%! assert_equal (mdl1.SSE, 11.6666666666667, -1e-9);
%! assert_equal (t.SumSq(1), 2564.33333333333, -1e-8);
%! assert_equal (t.DF(1), 3);
%! assert_equal (t.F(1), 512.866666666667, -1e-7);
%! assert_equal (t.pValue(1), 1.45297790506045e-08, -1e-8);
%! assert_equal (t.SumSq(2), 11.6666666666667, -1e-9);
%! assert_equal (t.DF(2), 7);
%! assert_equal (t.MeanSq(2), 1.66666666666667, -1e-9);

%!test
%! ## duplicate group column, so type 3 should show zero DF
%! G1 = categorical ([1;1;1;2;2;2;3;3;3;3]);
%! G2 = G1;
%! y1  = [10;12;11;20;22;19;15;17;16;14];
%! T = table (G1, G2, y1);
%! mdl1 = fitlm (T, 'y1 ~ G1 + G2');
%! t = anova (mdl1, 'components', 3);
%! assert_equal (t.Properties.RowNames, {'G1'; 'G2'; 'Error'});
%! assert_equal (t.SumSq(1), 0, -1e-9);
%! assert_equal (t.DF(1), 0);
%! assert (isnan (t.F(1)));
%! assert (isnan (t.pValue(1)));
%! assert_equal (t.SumSq(2), 0, -1e-9);
%! assert_equal (t.DF(2), 0);
%! assert_equal (t.SumSq(3), 11.6666666666667, -1e-9);
%! assert_equal (t.DF(3), 7);
%! assert_equal (t.MeanSq(3), 1.66666666666667, -1e-9);

%!test
%! ## just an intercept, so only the error row shows up
%! y1 = [5.1;4.9;5.2;4.8;5.1;4.9;5.2;4.8;5.1;4.9];
%! T = table (y1);
%! mdl1 = fitlm (T, 'y1 ~ 1');
%! t = anova (mdl1, 'components', 3);
%! assert_equal (t.Properties.RowNames, {'Error'});
%! assert_equal (t.SumSq(1), 0.22, -1e-9);
%! assert_equal (t.DF(1), 9);
%! assert_equal (t.MeanSq(1), 0.0244444444444444, -1e-9);

%!test
%! ## type 3 with a robust fit
%! G = categorical ([1;1;1;1;2;2;2;2;3;3;3;3]);
%! x = (1:12)' / 2;
%! y1 = 5 + 2*(G=='2') + 4*(G=='3') + 0.3*x;
%! y1(10) = y1(10) + 20;
%! T = table (G, x, y1);
%! mdl1 = fitlm (T, 'y1 ~ G + x', 'RobustOpts', 'on');
%! t = anova (mdl1, 'components', 3);
%! assert_equal (t.SumSq(1), 3.35664335664336, -1e-8);
%! assert_equal (t.F(1), 0.0801017164653529, -1e-8);
%! assert_equal (t.pValue(1), 0.923753304031669, -1e-9);
%! assert_equal (t.SumSq(2), 0.337499999999999, -1e-8);
%! assert_equal (t.F(2), 0.0161079545454545, -1e-8);
%! assert_equal (t.pValue(2), 0.902138027698746, -1e-9);
%! assert_equal (t.SumSq(3), 167.619047619048, -1e-7);
%! assert_equal (t.DF(3), 8);
%! assert_equal (t.MeanSq(3), 20.9523809523809, -1e-8);

%!test
%! ## repeated x values, so lack of fit rows should show up
%! x = [1;1;1;2;3;3;4;5;5;5;6;7];
%! y1 = [10.1;9.9;10.3;12.0;15.2;14.8;17.5;20.1;19.9;20.3;22.0;25.5];
%! mdl1 = fitlm (x, y1);
%! t = anova (mdl1, 'summary');
%! assert_equal (t.Properties.RowNames, ...
%!   {'Total'; 'Model'; 'Residual'; '. Lack of fit'; '. Pure error'});
%! assert_equal (t.SumSq(3), 1.03026642984015, -1e-9);
%! assert_equal (t.DF(3), 10);
%! assert_equal (t.SumSq(4), 0.790266429840146, -1e-9);
%! assert_equal (t.DF(4), 5);
%! assert_equal (t.F(4), 3.2927767910006, -1e-9);
%! assert_equal (t.pValue(4), 0.108434644423991, -1e-9);
%! assert_equal (t.SumSq(5), 0.24, -1e-9);
%! assert_equal (t.DF(5), 5);
%! assert_equal (t.MeanSq(5), 0.0480000000000001, -1e-9);

%!test
%! ## same as above but weighted
%! x = repmat ((1:4)', 6, 1);
%! G = categorical (repmat ([1;1;2;2], 6, 1));
%! y1 = 3 + 0.7*x + 2*(G=='2') + ...
%!     [0.1;-0.1;0.2;-0.2;0.15;-0.15;0.1;-0.1;-0.05;0.05;0.2;-0.2; ...
%!      0.1;-0.1;-0.15;0.15;0.2;-0.2;0.05;-0.05;-0.1;0.1;0.15;-0.15];
%! w = 1 + mod ((0:23)', 3);
%! T = table (x, G, y1);
%! mdl1 = fitlm (T, 'y1 ~ x + G', 'Weights', w);
%! t = anova (mdl1, 'summary');
%! assert_equal (t.Properties.RowNames, ...
%!   {'Total'; 'Model'; 'Residual'; '. Lack of fit'; '. Pure error'});
%! assert_equal (t.SumSq(3), 0.626927083333335, -1e-9);
%! assert_equal (t.DF(3), 21);
%! assert_equal (t.SumSq(4), 0.00880208333333332, -1e-9);
%! assert_equal (t.DF(4), 1);
%! assert_equal (t.F(4), 0.284799460734748, -1e-9);
%! assert_equal (t.pValue(4), 0.599454252702211, -1e-9);
%! assert_equal (t.SumSq(5), 0.618125000000001, -1e-9);
%! assert_equal (t.DF(5), 20);
%! assert_equal (t.MeanSq(5), 0.0309062500000001, -1e-9);

%!test
%! w = (1:n)';
%! mdl0 = fitlm (X, y, 'Weights', w);
%! m = step (mdl0, 'Upper', 'quadratic', 'Verbose', 0);
%! assert_equal (m.NumObservations, 20);
%! assert_equal (m.NumCoefficients, 3);
%! assert_equal (m.DFE, 17);
%! assert_equal (m.SSE, 4.16523310465, 1e-8);
%! assert_equal (m.CoefficientNames, {'(Intercept)','x1','x2'});
%! assert_equal (m.Coefficients.Estimate, [0.080321; 2.6483; -0.98452], 1e-4);
%! assert_equal (m.Formula.LinearPredictor, '1 + x1 + x2');

%!test
%! mdl0 = fitlm (X, y, 'Exclude', [3 7 15]);
%! m = step (mdl0, 'Upper', 'quadratic', 'Verbose', 0);
%! assert_equal (m.NumObservations, 17);
%! assert_equal (m.NumCoefficients, 3);
%! assert_equal (m.DFE, 14);
%! assert_equal (m.SSE, 0.340968585251, 1e-8);
%! assert_equal (m.CoefficientNames, {'(Intercept)','x1','x2'});
%! assert_equal (m.Coefficients.Estimate, [0.13263; 2.3566; -0.97211], 1e-4);
%! assert_equal (m.Formula.LinearPredictor, '1 + x1 + x2');

%!test
%! mdl0 = fitlm (X, y, 'Intercept', false);
%! m = step (mdl0, 'Upper', 'quadratic', 'Verbose', 0);
%! assert_equal (m.NumObservations, 20);
%! assert_equal (m.NumCoefficients, 2);
%! assert_equal (m.DFE, 18);
%! assert_equal (m.SSE, 0.410934843408, 1e-8);
%! assert_equal (m.Formula.HasIntercept, false);
%! assert_equal (m.CoefficientNames, {'x1','x2'});
%! assert_equal (m.Coefficients.Estimate, [2.9614; -0.99725], 1e-4);
%! assert_equal (m.Formula.LinearPredictor, 'x1 + x2');

%!test
%! load hald
%! Xh = ingredients;
%! yh = heat;
%! mdlh = fitlm (Xh, yh);
%! assert_equal (mdlh.Coefficients.Estimate, ...
%!   [62.405369299918; 1.55110264750845; 0.510167579684912; ...
%!    0.101909403579662; -0.144061029071018], 1e-9);
%! assert_equal (mdlh.Coefficients.SE, ...
%!   [70.0709592085362; 0.744769867130993; 0.72378800183518; ...
%!    0.754709045051309; 0.70905206344651], 1e-9);
%! assert_equal (mdlh.Coefficients.tStat, ...
%!   [0.890602469336764; 2.0826603169159; 0.704857746178952; ...
%!    0.135031379639465; -0.203174120065001], 1e-8);
%! assert_equal (mdlh.Coefficients.pValue, ...
%!   [0.399133563385561; 0.0708216874297252; 0.500901103474289; ...
%!    0.895922690510107; 0.844071473291884], 1e-8);
%! assert_equal (mdlh.CoefficientNames, {'(Intercept)', 'x1', 'x2', 'x3', 'x4'});
%! assert_equal (mdlh.NumCoefficients,          5);
%! assert_equal (mdlh.NumEstimatedCoefficients, 5);
%! assert_equal (mdlh.DFE,                      8);
%! assert_equal (mdlh.SSE, 47.863639350499,   1e-8);
%! assert_equal (mdlh.SSR, 2667.89943757258,  1e-6);
%! assert_equal (mdlh.SST, 2715.76307692308,  1e-6);
%! assert_equal (mdlh.MSE,  5.98295491881254, 1e-9);
%! assert_equal (mdlh.RMSE, 2.44600795559061, 1e-9);
%! assert_equal (mdlh.Rsquared.Ordinary, 0.98237562040768, 1e-9);
%! assert_equal (mdlh.Rsquared.Adjusted, 0.97356343061152, 1e-9);
%! assert_equal (mdlh.LogLikelihood, -26.918344895826, 1e-8);
%! assert_equal (mdlh.ModelCriterion.AIC,  63.8366897916521, 1e-7);
%! assert_equal (mdlh.ModelCriterion.AICc, 72.4081183630806, 1e-7);
%! assert_equal (mdlh.ModelCriterion.BIC,  66.6614365789598, 1e-7);
%! assert_equal (mdlh.ModelCriterion.CAIC, 71.6614365789598, 1e-7);
%! assert_equal (mdlh.Fitted(1:5), ...
%!   [78.4952395815018; 72.7887993002909; 105.970937532083; ...
%!    89.3271002550427; 95.649244438227], 1e-8);
%! assert_equal (mdlh.Residuals.Raw(1:5), ...
%!   [0.00476041849822195; 1.51120069970906; -1.67093753208295; ...
%!    -1.72710025504266; 0.250755561773033], 1e-8);
%! assert_equal (mdlh.Residuals.Pearson(1:5), ...
%!   [0.00194619910672879; 0.617823297040002; -0.683128412670876; ...
%!    -0.706089385807266; 0.102516249466771], 1e-8);
%! assert_equal (mdlh.Residuals.Studentized(1:5), ...
%!   [0.00271470565323249; 0.734526653667679; -1.05809320265782; ...
%!    -0.824036396702643; 0.119767490249399], 1e-8);
%! assert_equal (mdlh.Residuals.Standardized(1:5), ...
%!   [0.00290214088954622; 0.756624558354514; -1.05027405557414; ...
%!    -0.841081414787206; 0.127905848829164], 1e-8);
%! assert_equal (mdlh.Diagnostics.Leverage(1:5), ...
%!   [0.550284813713987; 0.333242829857405; 0.576942476415795; ...
%!    0.29523667959374; 0.357601364034465], 1e-8);
%! assert_equal (mdlh.Diagnostics.CooksDistance(1:5), ...
%!   [2.06118491039641e-06; 0.0572247602223712; 0.300862709270433; ...
%!    0.0592697490074745; 0.00182140011900327], 1e-8);
%! assert_equal (mdlh.Diagnostics.Dffits(1:5), ...
%!   [0.00300294746488125; 0.519283016757055; -1.23563576459509; ...
%!    -0.533347058289184; 0.0893585735072963], 1e-7);
%! assert_equal (mdlh.Diagnostics.S2_i(1:5), ...
%!   [6.83765556564705; 6.34835899957971; 5.89485540180634; ...
%!    6.23302709557492; 6.82367982420565], 1e-7);
%! assert_equal (mdlh.Diagnostics.CovRatio(1:5), ...
%!   [4.33530738335252; 2.01725612858557; 2.19476339013102; ...
%!    1.74129811023362; 3.00406926806094], 1e-6);
%! assert_equal (coefCI (mdlh), ...
%!   [-99.178552392689 223.989290992525; -0.166339745871082 3.26854504088797; ...
%!    -1.15889054555817 2.179225704928; -1.63845277518465 1.84227158234397; ...
%!    -1.77913801945372 1.49101596131168], 1e-7);
%! assert_equal (coefCI (mdlh, 0.1), ...
%!   [-67.8949453842232 192.705683984059; 0.166167302672858 2.93603799234403; ...
%!    -0.835750978716108 1.85608613808593; -1.30150832005232 1.50532712721164; ...
%!    -1.46257740216021 1.17445534401817], 1e-7);
%! [p, F, r] = coefTest (mdlh);
%! assert_equal (p, 4.75618174559791e-07, 1e-12);
%! assert_equal (F, 111.479171821258, 1e-6);
%! assert_equal (r, 4);
%! [dw, pdw] = dwtest (mdlh);
%! assert_equal (dw,  0.842123108585363, 1e-9);
%! assert_equal (pdw, 2.05259693286049,  1e-8);

%!test
%! load hald
%! Xq = ingredients;
%! yq = heat;
%! mdlq = fitlm (Xq, yq, 'purequadratic');
%! assert_equal (mdlq.Coefficients.Estimate, ...
%!   [-210.864812527187; 4.01775196981369; 5.27927179495849; 4.98703005469684; ...
%!    1.18967556414545; -0.00475259718542981; -0.0278555063988026; ...
%!    -0.0885225308739459; 0.0125632231921875], 1e-8);
%! assert_equal (mdlq.Coefficients.SE, ...
%!   [62.2492955088454; 0.709968572684776; 0.928573133899706; 1.15325562045609; ...
%!    0.559753134550446; 0.0119258668009245; 0.00570193167405062; ...
%!    0.0224063172978396; 0.0040139580330777], 1e-8);
%! assert_equal (mdlq.Coefficients.tStat, ...
%!   [-3.38742488253901; 5.6590560827508; 5.68535918413584; 4.3243058748108; ...
%!    2.12535757410437; -0.398511677579811; -4.88527537528597; ...
%!    -3.95078449069724; 3.12988404180068], 1e-7);
%! assert_equal (mdlq.Coefficients.pValue, ...
%!   [0.0275955163014823; 0.00480588842249351; 0.00472567672152044; ...
%!    0.0124052274432915; 0.100732674864514; 0.710609610173011; ...
%!    0.00812965117162355; 0.0168069450605605; 0.0351891921636521], 1e-7);
%! assert_equal (mdlq.CoefficientNames, ...
%!   {'(Intercept)', 'x1', 'x2', 'x3', 'x4', 'x1^2', 'x2^2', 'x3^2', 'x4^2'});
%! assert_equal (mdlq.NumCoefficients, 9);
%! assert_equal (mdlq.DFE, 4);
%! assert_equal (mdlq.Rsquared.Ordinary, 0.998060528051984, 1e-9);
%! assert_equal (mdlq.Rsquared.Adjusted, 0.994181584155951, 1e-9);
%! assert_equal (mdlq.SSE, 5.26714630515049, 1e-7);
%! assert_equal (mdlq.SSR, 2710.49593061793, 1e-6);
%! assert_equal (mdlq.SST, 2715.76307692308, 1e-6);
%! Xnewq = mean (Xq, 1);
%! [ypredq, yciq] = predict (mdlq, Xnewq);
%! assert_equal (ypredq, 101.90428629036, 1e-7);
%! assert_equal (yciq, [97.996264336243, 105.812308244477], 1e-6);
%! [ypredq2, yciq2] = predict (mdlq, Xnewq, 'Alpha', 0.01);
%! assert_equal (ypredq2, 101.90428629036, 1e-7);
%! assert_equal (yciq2, [95.4237317847484, 108.384840795971], 1e-6);
%! yfeq = feval (mdlq, Xnewq(1), Xnewq(2), Xnewq(3), Xnewq(4));
%! assert_equal (yfeq, 101.90428629036, 1e-7);
%! ysimq = random (mdlq, Xnewq);
%! assert_equal (isscalar (ysimq), true);
%! assert_equal (isnumeric (ysimq), true);
%! mdlq2 = removeTerms (mdlq, 'x1^2');
%! assert_equal (mdlq2.Coefficients.Estimate, ...
%!   [180.815670525319; 0.228012666917658; -2.0654020024496; -1.90840414992816; ...
%!    -0.0108001957000414; 0.0257318558607856; 0.00699693273966028], 1e-8);
%! assert_equal (mdlq2.CoefficientNames, ...
%!   {'(Intercept)', 'x2', 'x3', 'x4', 'x2^2', 'x3^2', 'x4^2'});
%! assert_equal (mdlq2.NumCoefficients, 7);
%! mdlq3 = addTerms (mdlq2, 'x1^2');
%! assert_equal (mdlq3.Coefficients.Estimate, mdlq.Coefficients.Estimate, 1e-8);
%! assert_equal (mdlq3.CoefficientNames, mdlq.CoefficientNames);
%! assert_equal (mdlq3.SSE, 5.26714630515049, 1e-7);

%!test
%! load hald
%! Xr = ingredients;
%! yr = heat;
%! mdlr = fitlm (Xr, yr, 'RobustOpts', 'bisquare');
%! assert_equal (mdlr.Coefficients.Estimate, ...
%!   [60.0897358816096; 1.57529551556915; 0.532199192097796; ...
%!    0.133455378556458; -0.120521170556001], 1e-8);
%! assert_equal (mdlr.Coefficients.SE, ...
%!   [75.8175597390933; 0.805849306629754; 0.783146694256936; ...
%!    0.816603608044244; 0.767202244491812], 1e-8);
%! assert_equal (mdlr.Coefficients.tStat, ...
%!   [0.792556976093573; 1.95482642053437; 0.679565138945976; ...
%!    0.163427368238162; -0.157091785668371], 1e-7);
%! assert_equal (mdlr.Coefficients.pValue, ...
%!   [0.450897370203866; 0.0863457969332376; 0.515957116726031; ...
%!    0.874235088124976; 0.879064839096153], 1e-7);
%! assert_equal (mdlr.CoefficientNames, {'(Intercept)', 'x1', 'x2', 'x3', 'x4'});
%! assert_equal (is_function_handle (mdlr.Robust.RobustWgtFun), true);
%! assert_equal (mdlr.Robust.Tune, 4.685, 1e-10);
%! assert_equal (size (mdlr.Robust.Weights), [13, 1]);
%! assert_equal (isnumeric (mdlr.Robust.Weights), true);
%! assert_equal (mdlr.SSE, 56.0362670671825, 1e-6);
%! assert_equal (mdlr.MSE, 7.00453338339782, 1e-8);
%! assert_equal (mdlr.RMSE, 2.64660790133291, 1e-8);
%! assert_equal (mdlr.Rsquared.Ordinary, 0.97929734395902, 1e-9);
%! assert_equal (mdlr.Rsquared.Adjusted, 0.96894601593853, 1e-9);
%! assert_equal (mdlr.DFE, 8);
%! H = [0 1 -1 0 0];
%! [p1, F1, r1] = coefTest (mdlr, H);
%! assert_equal (p1, 0.00308748318894346, 1e-11);
%! assert_equal (F1, 17.4568343157849, 1e-7);
%! assert_equal (r1, 1);
%! [pd1, dw1] = dwtest (mdlr, 'exact', 'both');
%! assert_equal (pd1, 0.844119247360191, 1e-9);
%! assert_equal (dw1, 2.05387711905232, 1e-8);
%! [pd2, dw2] = dwtest (mdlr, 'approximate', 'right');
%! assert_equal (pd2, 0.425180546504485, 1e-9);
%! assert_equal (dw2, 2.05387711905232, 1e-8);
%! Xnewr = mean (Xr, 1);
%! [ypredr, ycir] = predict (mdlr, Xnewr, 'Simultaneous', true);
%! assert_equal (ypredr, 95.4263340097424, 1e-7);
%! assert_equal (ycir, [92.2744598719279, 98.5782081475569], 1e-6);
%! [ypredr2, ycir2] = predict (mdlr, Xnewr, 'Simultaneous', true, 'Alpha', 0.1);
%! assert_equal (ypredr2, 95.4263340097424, 1e-7);
%! assert_equal (ycir2, [92.7161333049321, 98.1365347145527], 1e-6);

%!warning <addTerms: There are no new terms among the terms you specified.> ...
%! m = addTerms (mdl, 'x1');

%!warning <removeTerms: No specified terms appear in the model.> ...
%! m = removeTerms (mdl, 'x1:x2');

%!error <Unknown option 'NotAKey'> fitlm (X, y, 'NotAKey', 1)
%!error <VarNames must have 3 elements> fitlm (X, y, 'VarNames', {'a','b','c','d'})
%!error <Terms matrix must have 2 or 3 columns> fitlm (X, y, [1 2 3 4; 5 6 7 8])
%!error <Last column of terms matrix must be all zeros> fitlm (X, y, [1 2 1; 0 1 1])
%!error <No observations remain> fitlm (NaN (5, 2), NaN (5, 1))
%!error <No observations remain> fitlm (NaN (3, 2), [1; 2; 3])
%!error <No observations remain> fitlm ([1 2; 3 4; 5 6], NaN (3, 1))
%!error <No observations remain> fitlm (X, y, 'Exclude', (1:n)')
%!error <Not enough input arguments> fitlm ()
%!error <Predictor variables must be numeric> fitlm ('hello', y)
%!error <Predictor variables must be numeric> fitlm ({'a';'b'}, [1; 2])
%!error <Y argument is required> fitlm (X)
%!error <Y argument is required> fitlm (X, 'Weights', [1;1;1])
%!error <Predictor and response variables must have the same length> fitlm (X, [1; 2])
%!error <Predictor and response variables must have the same length> fitlm (X, [1 2])
%!error <indexing is not supported> mdl(1)
%!error <indexing is not supported> mdl {1}
%!error <unknown option> predict (mdl, [0.5 0.25], 'BadOption', 1)
%!error <Alpha must be a scalar> predict (mdl, [0.5 0.25], 'Alpha', -0.1)
%!error <Alpha must be a scalar> predict (mdl, [0.5 0.25], 'Alpha', 1.5)
%!error <Alpha must be a scalar> predict (mdl, [0.5 0.25], 'Alpha', [0.01 0.05])
%!error <Prediction must be> predict (mdl, [0.5 0.25], 'Prediction', 'bad')
%!error <Xnew must have 2 columns> predict (mdl, ones (3, 5))
%!error <Xnew must have 2 columns> predict (mdl, ones (3, 1))
%!error <missing predictor> predict (mdl, table ([1;2], 'VariableNames', {'z'}))
%!error <Not enough input arguments> random (mdl)
%!error <Too many input arguments> random (mdl, [0.5, 0.25], 'extra')
%!error <Xnew must have 2 columns> random (mdl, ones (3, 5))
%!error <Xnew must have 2 columns> random (mdl, [])
%!error <Not enough input arguments> feval (mdl)
%!error <Incorrect number of input arguments> feval (mdl, [0.5; 1.0], [0.25; 1.0], [0.1; 0.2])
%!error <Predictor data matrix must have 2 columns> feval (mdl, ones (3, 1))
%!error <All input arguments must be the same size> feval (mdl, [0.5; 1.0; 0.2], [0.25; 1.0])
%!error <X does not contain one or more predictor> feval (mdl, table ([1; 2], 'VariableNames', {'z'}))
%!error <Predictor data matrix must have 2 columns> feval (mdl, [])
%!error <is not categorical> feval (mdl, '0.5', 0.25)
%!error <too many inputs> coefCI (mdl, 0.05, 'extra')
%!error <Value must be less than or equal to 1> coefCI (mdl, 1.5)
%!error <Value must be greater than or equal to 0> coefCI (mdl, -0.1)
%!error <Value must be a scalar> coefCI (mdl, [0.01 0.05])
%!error <Value must be greater than or equal to 0> coefCI (mdl, NaN)
%!error <Value must be a scalar> coefCI (mdl, 'abc')
%!error <H must be a 1-by-3 numeric matrix> coefTest (mdl, [1 0])
%!error <H must be a 1-by-3 numeric matrix> coefTest (mdl, 'abc')
%!error <C must be a numeric vector> coefTest (mdl, [0 1 0], 'abc')
%!error <H must be a 1-by-3 numeric matrix> coefTest (mdl, [0 1 0; 0 0 1], [1])
%!error <H is not full rank> coefTest (mdl, [0 NaN 0])
%!error <Too many input arguments> coefTest (mdl, [0 1 0], 0, 'extra')
%!error <too many outputs> [a, b, c, d] = coefTest (mdl)
%!error <The METHOD argument must be> dwtest (mdl, 'badmethod', 'both')
%!error <The METHOD argument must be> dwtest (mdl, 123, 'both')
%!error <Too many input arguments> dwtest (mdl, 'exact', 'both', 'extra')
%!error <too many outputs> [a, b, c] = dwtest (mdl)
%!error <Not enough input arguments> addTerms (mdl)
%!error <too many inputs> addTerms (mdl, 'x1:x2', 'extra')
%!error <Unrecognized variable> addTerms (mdl, 'z')
%!error <Unrecognized variable> addTerms (mdl, 'X1')
%!error <Unrecognized variable> addTerms (mdl, 'x1:z')
%!error <Unrecognized variable> addTerms (mdl, 'x1*z')
%!error <Terms matrix must have> addTerms (mdl, [1, 1, 1, 0])
%!error <Terms matrix must have> addTerms (mdl, [])
%!error <Model update specification> addTerms (mdl, {})
%!error <Not enough input arguments> removeTerms (mdl)
%!error <too many inputs> removeTerms (mdl, 'x1', 'extra')
%!error <Unrecognized variable> removeTerms (mdl, 'z1')
%!error <Terms matrix must have 3 columns> removeTerms (mdl, [0 1 0 0])
%!error <Terms matrix must have 3 columns> removeTerms (mdl, [])
%!error <Model update specification> removeTerms (mdl, {'x1'})
%!error <Bad residuals plot type> plotResiduals (mdl, 'badtype')
%!error <invalid ResidualType> plotResiduals (mdl, 'fitted', 'ResidualType', 'bad')
%!error <Bad diagnostics plot type> plotDiagnostics (mdl, 'badtype')
%!error <unrecognized property> plotDiagnostics (mdl, 'leverage', 'BadProp', 1)
%!error <Wrong number of arguments> plotEffects (mdl, 'extra')
%!error <Wrong number of arguments> plotEffects (mdl, 'a', 'b')
%!error <Model has no predictors> plotEffects (fitlm (X(:,1), y, 'constant'))
%!error <unrecognised RobustWgtFun> fitlm (X, y, 'RobustOpts', 'notarealfunction')
%!error <invalid RobustOpts value> fitlm (X, y, 'RobustOpts', 42)
%!error <Not enough input arguments> plotAdjustedResponse (mdl)
%!error <is not a variable for this fit> plotAdjustedResponse (mdl, 'z')
%!error <is the response in this model> plotAdjustedResponse (mdl, 3)
%!error <This model only contains> plotAdjustedResponse (mdl, 99)
%!error <Variable must be specified as a name or a positive integer> plotAdjustedResponse (mdl, 1.5)
%!error <unrecognized property> plotAdjustedResponse (mdl, 'x1', 'BadOption', 5)
%!error <Bad coefficient number> plotAdded (mdl, 99)
%!error <Bad coefficient name> plotAdded (mdl, 'NotACoef')
%!error <unrecognized property> plotAdded (mdl, 2, 'BadOpt', 5)
%!error <Bad coefficient number> mdl0 = fitlm (ones (n, 1), y, 'Intercept', false); plotAdded (mdl0)
%!error <Too many input arguments> plot (mdl, 'extra')
%!error <Not enough input arguments> plotInteraction (mdl)
%!error <Not enough input arguments> plotInteraction (mdl, 'x1')
%!error <PTYPE must be> plotInteraction (mdl, 'x1', 'x2', 'badtype')
%!error <Too many input arguments> plotInteraction (mdl, 'x1', 'x2', 'effects', 'extra')
%!error <is not a variable for this fit> plotInteraction (mdl, 'z', 'x2')
%!error <is not a variable for this fit> plotInteraction (mdl, 'x1', 'z')
%!error <This model only contains> plotInteraction (mdl, 99, 'x2')
%!error <Variable must be specified as a name or a positive integer> plotInteraction (mdl, 1.5, 'x2')
%!error <is the response in this model> plotInteraction (mdl, 'y', 'x2')
%!error <is the response in this model> plotInteraction (mdl, 'x1', 'y')
%!error <VAR1 and VAR2 must be different variables> plotInteraction (mdl, 'x1', 'x1')
%!error <too many inputs> compact (mdl, 'extra')
%!error <anova: too many input arguments.> anova (mdl, 'components', 'h', 'extra')
%!error <anova: ANOVATYPE must be 'summary' or 'components'.> anova (mdl, 'bogus')
%!error <anova: SSTYPE can only be specified with ANOVATYPE 'components'.> anova (mdl, 'summary', 2)
%!error <anova: SSTYPE must be 1, 2, or 3.> anova (mdl, 'components', 4)
%!error <The STEP method is not available with a robust fit> ...
%! mdl0 = fitlm (X, y, 'RobustOpts', 'on'); step (mdl0)
%!error <Name-Value arguments must be in pairs> step (mdl, 'Verbose')

## A factor reaching the model only through an interaction or a power has a
## column of the terms matrix without ever being a coefficient, so anything
## indexed by those columns must be built from their own names.
%!shared tf, tv, tr, ttbl
%! tf = [1;2;3;4;5;6;7;8;9;10;11;12];
%! tv = [2;1;4;3;6;5;8;7;10;9;12;11];
%! tr = [3.1;4.2;5.3;6.4;7.5;8.6;9.7;10.8;11.9;13.0;14.1;15.2];
%! ttbl = table (tf, tv, tr, 'VariableNames', {'u', 'v', 'resp'});

%!test  # a power keeps its exponent in the variable's own column
%! m = fitlm (ttbl, 'resp ~ 1 + u*v + u^2 + v^2');
%! assert_equal (char (m.Formula), 'resp ~ 1 + u*v + u^2 + v^2');
%! assert_equal (m.Formula.Terms, ...
%!               [0 0 0; 1 0 0; 0 1 0; 1 1 0; 2 0 0; 0 2 0]);
%! assert_equal (m.Formula.TermNames, ...
%!               {'(Intercept)'; 'u'; 'v'; 'u:v'; 'u^2'; 'v^2'});

%!test  # an interaction survives when one factor is not a main effect
%! m = fitlm (ttbl, 'resp ~ 1 + u + u:v');
%! assert_equal (char (m.Formula), 'resp ~ 1 + u + u:v');
%! assert_equal (m.Formula.Terms, [0 0 0; 1 0 0; 1 1 0]);
%! assert_equal (m.Formula.TermNames, {'(Intercept)'; 'u'; 'u:v'});
%! assert_equal (m.CoefficientNames, {'(Intercept)', 'u', 'u:v'});

%!test  # a power alone still names the variable it belongs to
%! m = fitlm (ttbl, 'resp ~ 1 + v^2 + u^2');
%! assert_equal (char (m.Formula), 'resp ~ 1 + u + v + u^2 + v^2');
%! assert_equal (m.Formula.NTerms, 5);
%! assert_equal (m.Formula.PredictorNames, {'u', 'v'});

## Dropping the intercept gives the first categorical predictor every one of
## its levels -- the cell-means parameterisation -- on both the formula and the
## keyword path.  MATLAB drops the reference level either way and so cannot fit
## the reference group; see the note in the class documentation.
%!shared cf, ch, cg, ct1, ct2
%! cf = [1;2;3;4;5;6;7;8;9;10;11;12];
%! ch = {'b';'c';'a';'b';'c';'a';'b';'c';'a';'b';'c';'a'};
%! cg = {'lo';'hi';'lo';'hi';'lo';'hi';'lo';'hi';'lo';'hi';'lo';'hi'};
%! ct1 = table (cf, ch, [3.1;4.2;5.3;6.4;7.5;8.6;9.7;10.8;11.9;13;14.1;15.2], ...
%!              'VariableNames', {'u', 'h2', 'resp'});
%! ct2 = table (ch, cg, [3.1;4.2;5.3;6.4;7.5;8.6;9.7;10.8;11.9;13;14.1;15.2], ...
%!              'VariableNames', {'h2', 'g', 'resp'});

%!test  # every level is coded, and the estimates are the group means
%! m = fitlm (ct1, 'resp ~ h2 - 1');
%! assert_equal (m.CoefficientNames, {'h2_b', 'h2_c', 'h2_a'});
%! assert_equal (m.Coefficients.Estimate', [8.05, 9.15, 10.25], 1e-12);
%! assert_equal (m.Rsquared.Ordinary > 0, true);

%!test  # the keyword path codes it the same way
%! m = fitlm (ct1, 'linear', 'Intercept', false);
%! assert_equal (m.CoefficientNames, {'u', 'h2_b', 'h2_c', 'h2_a'});
%! assert_equal (m.NumCoefficients, m.NumEstimatedCoefficients);

%!test  # a second categorical stays reference coded, keeping the design full rank
%! for spec = {'resp ~ h2 + g - 1', 'linear'}
%!   if (strcmp (spec{1}, 'linear'))
%!     m = fitlm (ct2, spec{1}, 'Intercept', false);
%!   else
%!     m = fitlm (ct2, spec{1});
%!   endif
%!   assert_equal (m.CoefficientNames, {'h2_b', 'h2_c', 'h2_a', 'g_hi'});
%!   assert_equal (m.NumCoefficients, m.NumEstimatedCoefficients);
%! endfor

%!test  # an intercept still takes the reference level with it
%! m = fitlm (ct1, 'resp ~ 1 + h2');
%! assert_equal (m.CoefficientNames, {'(Intercept)', 'h2_c', 'h2_a'});
%! m = fitlm (ct1, 'linear');
%! assert_equal (m.CoefficientNames, {'(Intercept)', 'u', 'h2_c', 'h2_a'});

## A predictor that groups its observations is coded the same way whatever its
## column type, and VariableInfo reports it as categorical with its levels in
## the order the design codes them.
%!shared vu, vg, vb, vk, vr
%! vu = [1;2;3;4;5;6;7;8;9;10;11;12];
%! vg = {'lo';'hi';'lo';'hi';'lo';'hi';'lo';'hi';'lo';'hi';'lo';'hi'};
%! vb = logical ([0;1;0;1;0;1;0;1;0;1;0;1]);
%! vk = [5;3;5;3;5;3;5;3;5;3;5;3];
%! vr = [3.1;4.2;5.3;6.4;7.5;8.6;9.7;10.8;11.9;13;14.1;15.2];

%!test  # a logical column is coded by value, with false as the reference
%! t = table (vu, vb, vr, 'VariableNames', {'u', 'b', 'resp'});
%! m = fitlm (t, 'resp ~ 1 + u + b');
%! assert_equal (m.CoefficientNames, {'(Intercept)', 'u', 'b_1'});
%! assert_equal (m.VariableInfo.IsCategorical, [false; true; false]);
%! assert_equal (m.VariableInfo.Class, {'double'; 'logical'; 'double'});
%! assert_equal (m.VariableInfo.Range{2}, [false, true]);
%! ## and the keyword path codes it identically
%! assert_equal (fitlm (t, 'linear').CoefficientNames, ...
%!               {'(Intercept)', 'u', 'b_1'});

%!test  # 'CategoricalVars' is honoured on the formula path too
%! t = table (vu, vk, vr, 'VariableNames', {'u', 'k', 'resp'});
%! m = fitlm (t, 'resp ~ 1 + u + k', 'CategoricalVars', {'k'});
%! assert_equal (m.CoefficientNames, {'(Intercept)', 'u', 'k_5'});
%! assert_equal (m.VariableInfo.IsCategorical, [false; true; false]);
%! assert_equal (m.VariableInfo.Range{2}, [3, 5]);
%! ## without the declaration it stays numeric
%! m = fitlm (t, 'resp ~ 1 + u + k');
%! assert_equal (m.CoefficientNames, {'(Intercept)', 'u', 'k'});

%!test  # a string column groups like a cellstr one
%! t = table (vu, string (vg), vr, 'VariableNames', {'u', 'g', 'resp'});
%! m = fitlm (t, 'resp ~ 1 + u + g');
%! assert_equal (m.CoefficientNames, {'(Intercept)', 'u', 'g_hi'});
%! assert_equal (m.VariableInfo.IsCategorical, [false; true; false]);
%! assert_equal (class (m.VariableInfo.Range{2}), 'string');

%!test  # Range keeps each variable's own type and the design's level order
%! t = table (vu, vg, vr, 'VariableNames', {'u', 'g', 'resp'});
%! m = fitlm (t, 'resp ~ 1 + u + g');
%! assert_equal (m.VariableInfo.Range{2}, {'lo', 'hi'});
%! assert_equal (m.VariableInfo.Range{1}, [1, 12]);
%! t = table (vu, categorical (vg), vr, 'VariableNames', {'u', 'g', 'resp'});
%! m = fitlm (t, 'resp ~ 1 + u + g');
%! assert_equal (class (m.VariableInfo.Range{2}), 'categorical');
%! assert_equal (cellstr (m.VariableInfo.Range{2}), {'hi', 'lo'});
