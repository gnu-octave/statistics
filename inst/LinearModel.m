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
  ## @multitable @columnfractions 0.22 0.02 0.76
  ## @headitem Group @tab @tab Properties
  ##
  ## @item Coefficient estimates @tab @tab @code{Coefficients} (a table of
  ## estimates, standard errors, t-statistics, and p-values for each term),
  ## @code{CoefficientCovariance}, @code{CoefficientNames}, and the
  ## coefficient counts @code{NumCoefficients} and
  ## @code{NumEstimatedCoefficients}.
  ##
  ## @item Summary statistics of the fit @tab @tab @code{DFE},
  ## @code{Fitted}, @code{Residuals} (raw, Pearson, Studentized, and
  ## standardized), @code{Diagnostics} (leverage, Cook's distance, and other
  ## per-observation influence measures), @code{MSE}, @code{RMSE},
  ## @code{Rsquared} (ordinary and adjusted), @code{SSE}, @code{SSR},
  ## @code{SST}, @code{LogLikelihood}, @code{ModelCriterion} (AIC, BIC, etc.),
  ## and @code{ModelFitVsNullModel} (the F-test of the fitted model against an
  ## intercept-only model).
  ##
  ## @item Fitting method information @tab @tab @code{Robust}, which records
  ## the weighting function and tuning constant used when the model is fit by
  ## robust regression, and is empty for an ordinary least squares fit.
  ##
  ## @item Input data properties @tab @tab @code{Formula},
  ## @code{NumObservations}, @code{NumPredictors}, @code{NumVariables},
  ## @code{ObservationInfo} (which observations were used, excluded, missing,
  ## or weighted), @code{ObservationNames}, @code{PredictorNames},
  ## @code{ResponseName}, @code{VariableInfo}, @code{VariableNames}, and
  ## @code{Variables}.
  ## @end multitable
  ##
  ## A @code{LinearModel} object supports categorical predictors, which are
  ## automatically encoded internally as indicator (dummy) variables,
  ## observation weights for a weighted least squares fit, excluding specific
  ## observations from the fit, and robust regression using iteratively
  ## reweighted least squares.  Once fitted, the following methods are
  ## available on a @code{LinearModel} object:
  ##
  ## @multitable @columnfractions 0.2 0.02 0.78
  ## @headitem Method @tab @tab Description
  ##
  ## @item @code{predict} @tab @tab Predict responses at new predictor values
  ## given in a matrix or table, or reproduce the training fitted values when
  ## called with no new data.  Can also return pointwise or simultaneous
  ## confidence or prediction intervals alongside the point predictions.
  ##
  ## @item @code{feval} @tab @tab Predict responses given predictors as
  ## separate scalar or vector arguments (one per predictor variable) instead
  ## of a single matrix, so a @code{LinearModel} object can be evaluated the
  ## same way as a plain function handle.  Returns point predictions only.
  ##
  ## @item @code{random} @tab @tab Simulate new response values at new
  ## predictor locations by adding independent Gaussian noise, drawn from the
  ## estimated error variance @code{MSE}, to the fitted response.
  ##
  ## @item @code{coefCI} @tab @tab Return Wald confidence intervals for every
  ## fitted coefficient at a chosen significance level (default @math{0.05}).
  ##
  ## @item @code{coefTest} @tab @tab Test a linear hypothesis on the fitted
  ## coefficients.  With no arguments, tests the overall model F-test that
  ## all non-intercept coefficients are zero; a custom hypothesis can be
  ## given as a contrast matrix and, if needed, right-hand-side values.
  ## Returns the p-value, and optionally the F-statistic and its numerator
  ## degrees of freedom.
  ##
  ## @item @code{dwtest} @tab @tab Durbin-Watson test for first-order
  ## autocorrelation among the model residuals, with a choice of exact or
  ## approximate p-value computation and a one- or two-sided alternative.
  ##
  ## @item @code{addTerms} @tab @tab Return a new, refitted @code{LinearModel}
  ## with terms added to the current model specification, given as a
  ## Wilkinson formula fragment or a terms matrix.  Weights, excluded rows,
  ## and categorical encodings carry over automatically; the original model
  ## object is left unmodified.
  ##
  ## @item @code{removeTerms} @tab @tab Return a new, refitted
  ## @code{LinearModel} with terms removed from the current model
  ## specification, given as a Wilkinson formula fragment or a terms matrix.
  ## Weights, excluded rows, and categorical encodings carry over
  ## automatically; the original model object is left unmodified.
  ##
  ## @item @code{plotResiduals} @tab @tab Plot the model residuals.  Default
  ## is a probability density histogram; other supported plot types are
  ## @qcode{"fitted"}, @qcode{"caseorder"}, @qcode{"lagged"},
  ## @qcode{"probability"}, and @qcode{"observed"}.
  ##
  ## @item @code{plotDiagnostics} @tab @tab Plot per-observation influence
  ## diagnostics.  Default is leverage by observation row number; other
  ## supported plot types are @qcode{"cookd"}, @qcode{"covratio"},
  ## @qcode{"dfbetas"}, @qcode{"dffits"}, @qcode{"s2_i"}, and
  ## @qcode{"contour"} (standardized residuals against leverage with Cook's
  ## distance contours).
  ##
  ## @item @code{plotEffects} @tab @tab Plot the estimated main effect and
  ## 95% confidence interval of each predictor, evaluated between its
  ## observed minimum and maximum with all other predictors held at their
  ## observed means.
  ##
  ## @item @code{plotAdjustedResponse} @tab @tab Plot the fitted response
  ## against a single predictor, with the other predictors averaged out by
  ## averaging the fitted values over the observations used in the fit.
  ##
  ## @item @code{plotAdded} @tab @tab Plot the incremental effect of one or
  ## more terms on the response, after removing the effects of all other
  ## terms, along with the fitted line and its 95% confidence bounds.
  ## @end multitable
  ##
  ## Create a @code{LinearModel} object by using the @code{fitlm} function or
  ## the class constructor directly.
  ##
  ## @seealso{fitlm}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

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
    ## A structure with seven fields:
    ## @itemize
    ## @item @code{Start} - formula string of the starting model
    ## @item @code{Lower} - formula string of the lower-bound model; terms
    ##   listed here cannot be removed
    ## @item @code{Upper} - formula string of the upper-bound model; the
    ##   model cannot grow beyond this
    ## @item @code{Criterion} - criterion used, e.g. @qcode{'sse'}
    ## @item @code{PEnter} - threshold for adding a term
    ## @item @code{PRemove} - threshold for removing a term
    ## @item @code{History} - table with one row per step and columns
    ##   @code{Action}, @code{TermName}, @code{Terms}, @code{DF},
    ##   @code{delDF}, @code{FStat}, @code{PValue}
    ## @end itemize
    ## This structure is empty unless the model was fit using stepwise
    ## regression.  This property is read-only.
    ##
    ## @end deftp
    Steps = [];


    ## Input data properties

    ## -*- texinfo -*-
    ## @deftp {LinearModel} {property} Formula
    ##
    ## Model formula information
    ##
    ## A structure representing the model formula with fields including
    ## @code{ResponseName}, @code{LinearPredictor}, @code{PredictorNames},
    ## @code{TermNames}, @code{HasIntercept}, @code{Terms} (the terms
    ## matrix), and @code{InModel}.  This property is read-only.
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

  properties (Access = private, Hidden)

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

    ## Terms matrix from modelspec or lm_parse_modelspec
    TermsMatrix = [];

    ## Categorical level info for re-encoding in predict
    CatLevelInfo = [];

    ## Predictor names after categorical dummy expansion
    EncPredictorNames = {};

    ## Encoded predictor matrix (Path B only), cached for refit
    EncodedPredMatrix = [];

    ## Parsed NV options stored for refit operations
    OrigOpts = [];

  endproperties

  methods (Hidden)

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
      if (isempty (this.Robust))
        fprintf ("\n  Linear regression model:\n");
      else
        fprintf ("\n  Linear regression model (robust fit):\n");
      endif
      if (! isempty (this.Formula) && isstruct (this.Formula) ...
          && isfield (this.Formula, 'LinearPredictor'))
        fprintf ("      %s ~ %s\n", this.ResponseName, ...
                 this.Formula.LinearPredictor);
      endif

      if (! isempty (this.Coefficients))
        fprintf ("\n  Coefficients:\n\n");
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
          error (["LinearModel: () indexing is not supported.  " ...
                  "Use dot notation to access properties."]);
        case '{}'
          error (["LinearModel: {} indexing is not supported.  " ...
                  "Use dot notation to access properties."]);
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

  methods (Access = public)

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
    ## predictor data @var{X}.  Unless removed via the @qcode{"Intercept"}
    ## option, the fitted model contains a constant (intercept) term and one
    ## linear term for every column of @var{X}.
    ##
    ## @itemize
    ## @item
    ## @var{X} is an @math{NxP} numeric or logical matrix of predictor data,
    ## where rows correspond to observations and columns correspond to
    ## variables.  By default, the predictors are named @qcode{"x1"},
    ## @qcode{"x2"}, @dots{}, @qcode{"xP"}.
    ## @item
    ## @var{y} is an @math{Nx1} numeric or logical vector of response values,
    ## and must have the same number of observations (rows) as @var{X}.  By
    ## default, the response is named @qcode{"y"}.
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
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Value} @tab @tab @var{Description}
    ##
    ## @item @qcode{"constant"} @tab @tab Model contains only an intercept
    ## term.
    ##
    ## @item @qcode{"linear"} @tab @tab Model contains an intercept and one
    ## term for each predictor variable.  This is the default when
    ## @var{modelspec} is not specified.
    ##
    ## @item @qcode{"interactions"} @tab @tab Model contains an intercept, all
    ## linear terms, and all pairwise products of distinct predictor
    ## variables (no squared terms).
    ##
    ## @item @qcode{"purequadratic"} @tab @tab Model contains an intercept,
    ## all linear terms, and all squared terms.
    ##
    ## @item @qcode{"quadratic"} @tab @tab Model contains an intercept, all
    ## linear terms, all pairwise products of distinct predictor variables,
    ## and all squared terms.
    ##
    ## @item @qcode{"full"} @tab @tab Model contains an intercept and all
    ## terms up to and including the full @math{P}-way interaction of the
    ## predictor variables.
    ##
    ## @item terms matrix @tab @tab A @math{TxP} or @math{Tx(P+1)} numeric
    ## matrix, where @math{T} is the number of terms and @math{P} is the
    ## number of predictor variables.  Each row represents one term, and the
    ## value in column @math{j} is the exponent to which predictor @math{j}
    ## is raised in that term; a row of all zeros represents the intercept.
    ## If a @math{Tx(P+1)} matrix is supplied, its last column (representing
    ## the response variable) must be all zeros.
    ##
    ## @item Wilkinson formula @tab @tab A character vector of the form
    ## @qcode{"y ~ terms"} describing the response and predictor terms using
    ## Wilkinson notation.  For table input, the variable to the left of
    ## @qcode{"~"} is used as the response, overriding @var{resp_input}.
    ## @end multitable
    ##
    ## @code{@var{mdl} = LinearModel (@dots{}, @var{Name}, @var{Value},
    ## @dots{})} specifies additional options using one or more
    ## @qcode{Name-Value} pair arguments as described below.
    ##
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"Intercept"} @tab @tab A logical scalar indicating
    ## whether to include a constant (intercept) term in the model.  Default
    ## is @qcode{true}.  Ignored when @var{modelspec} is a Wilkinson formula.
    ##
    ## @item @qcode{"Weights"} @tab @tab A numeric vector of nonnegative
    ## observation weights, with one element per observation, used to fit a
    ## weighted least squares model.  Default is a vector of ones.
    ##
    ## @item @qcode{"Exclude"} @tab @tab A numeric or logical vector
    ## specifying observations to exclude from the fit, given as row indices
    ## or a logical mask.  Excluded observations, together with any
    ## observation containing a missing value, are recorded in
    ## @code{ObservationInfo} but do not contribute to the fit.
    ##
    ## @item @qcode{"CategoricalVars"} @tab @tab Specifies which predictor
    ## variables are treated as categorical, given as a vector of column
    ## indices, a logical vector, or a cell array of variable names.  Each
    ## categorical predictor with @math{L} categories is expanded into
    ## @math{L-1} indicator (dummy) variables, using the first category as
    ## the reference level.
    ##
    ## @item @qcode{"VarNames"} @tab @tab A cell array of character vectors
    ## naming the predictor and response variables, in order, with the
    ## response variable name last.  Only applies to matrix input, since
    ## table variables already carry their own names.
    ##
    ## @item @qcode{"ResponseVar"} @tab @tab A character vector naming the
    ## response variable, used to override the response variable name that
    ## would otherwise be used.
    ##
    ## @item @qcode{"PredictorVars"} @tab @tab A cell array of character
    ## vectors naming which variables in @var{tbl} to use as predictors.  By
    ## default, all variables other than the response variable are used.
    ##
    ## @item @qcode{"RobustOpts"} @tab @tab Selects ordinary least squares or
    ## robust regression fitting.  This value can be @qcode{"off"} (default,
    ## ordinary least squares), @qcode{"on"} (robust fitting using the
    ## @qcode{"bisquare"} weighting function), the name of one of the
    ## weighting functions below, a function handle for a custom weighting
    ## function, or a scalar structure with fields @qcode{RobustWgtFun} and
    ## @qcode{Tune} specifying the weighting function and its tuning
    ## constant.  Robust fitting uses Iteratively Reweighted Least Squares
    ## (IRLS), refitting the model with updated observation weights until the
    ## coefficients converge.  Supported weighting function names:
    ## @qcode{"andrews"}, @qcode{"bisquare"}, @qcode{"cauchy"},
    ## @qcode{"fair"}, @qcode{"huber"}, @qcode{"logistic"}, @qcode{"ols"},
    ## @qcode{"talwar"}, @qcode{"welsch"}, each with its own default tuning
    ## constant.
    ## @end multitable
    ##
    ## @var{mdl} is returned as a @code{LinearModel} object.  If
    ## @qcode{"RobustOpts"} is anything other than @qcode{"off"}, the returned
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
            error ('LinearModel: VarNames must have %d elements.', p_raw + 1);
          endif
          pred_names_raw = opts.VarNames(1:p_raw);
          resp_name      = opts.VarNames{end};
        else
          pred_names_raw = arrayfun (@(k) sprintf ('x%d', k), 1:p_raw, ...
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

        if (ischar (resp_input) && ! isempty (resp_input))
          resp_name = resp_input;
        elseif (isstring (resp_input) && ! isempty (resp_input))
          resp_name = char (resp_input);
        elseif (isnumeric (resp_input) && ! isempty (resp_input))
          resp_name = 'y';
          if (! isempty (opts.ResponseVar))
            resp_name = opts.ResponseVar;
          endif
          y_ext = double (resp_input(:));
        elseif (is_formula)
          tparts    = strsplit (modelspec, '~');
          resp_name = strtrim (tparts{1});
        else
          resp_name = col_names{end};
        endif

        if (! isempty (opts.PredictorVars))
          pred_names_raw = opts.PredictorVars;
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
          if (iscell (col) || isa (col, 'categorical'))
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
        error ('LinearModel: No observations remain after removing missing/excluded rows.');
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
          modelspec, 'model_matrix', tbl_sub);

        coef_names    = coef_names_raw(:)';
        y_sub         = y_full(subset_mask);
        n_coef        = size (X_design_sub, 2);
        has_intercept = any (strcmp (coef_names, '(Intercept)'));
        enc_names     = coef_names(! strcmp (coef_names, '(Intercept)'));

        cat_info.names  = {};
        cat_info.levels = {};
        if (istable (data))
          for j = 1:p_raw
            if (cat_logical(j))
              col = tbl_sub.(pred_names_raw{j});
              if (iscell (col))
                levels_j = unique (col);
              elseif (isa (col, 'categorical'))
                levels_j = categories (col);
              else
                levels_j = {};
              endif
              cat_info.names{end+1}  = pred_names_raw{j};
              cat_info.levels{end+1} = levels_j;
            endif
          endfor
        endif

        atomic_names = {};
        for t = 1:n_coef
          if (strcmp (coef_names{t}, '(Intercept)'))
            continue;
          endif
          factors_t = strsplit (coef_names{t}, ':');
          for f = 1:numel (factors_t)
            if (! any (strcmp (atomic_names, factors_t{f})))
              atomic_names{end+1} = factors_t{f};
            endif
          endfor
        endfor

        terms = zeros (n_coef, numel (atomic_names) + 1);
        for t = 1:n_coef
          if (strcmp (coef_names{t}, '(Intercept)'))
            continue;
          endif
          factors_t = strsplit (coef_names{t}, ':');
          for f = 1:numel (factors_t)
            col_t = find (strcmp (atomic_names, factors_t{f}));
            terms(t, col_t) = 1;
          endfor
        endfor

      else

        ## PATH B: Keyword / numeric terms matrix
        X_num_full     = zeros (n_total, p_raw);
        cat_str_levels = cell (1, p_raw);

        for j = 1:p_raw
          if (istable (data))
            col = tbl.(pred_names_raw{j});
            if (iscell (col))
              [cat_str_levels{j}, ~, ic] = unique (col);
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
              cat_str_levels{j} = cellstr (num2str (uvals(:)));
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

        [X_enc_sub, enc_names, cat_info] = lm_encode_categorical ( ...
          X_num_sub, cat_logical, pred_names_raw, cat_str_levels);
        p_enc  = size (X_enc_sub, 2);

        [terms, has_intercept, coef_names] = lm_parse_modelspec ( ...
          modelspec, enc_names, p_enc, opts.Intercept);
        n_coef = rows (terms);

        X_design_sub = LinearModel.lm_build_design (terms, X_enc_sub);

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
      CD_full  = NaN   (n_total, 1);
      Dff_full = NaN   (n_total, 1);
      S2i_full = NaN   (n_total, 1);
      CR_full  = NaN   (n_total, 1);
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
        for L = 2:numel (levels_c)
          dummy_names{end+1} = [base_nm, '_', char(levels_c{L})];
          dummy_bases{end+1} = base_nm;
        endfor
      endfor

      if (has_intercept)
        non_int = coef_names(! strcmp (coef_names, '(Intercept)'));
      else
        non_int = coef_names;
      endif

      disp_terms = {};
      for t = 1:numel (non_int)
        factors_t = strsplit (non_int{t}, ':');
        for f = 1:numel (factors_t)
          idx = find (strcmp (dummy_names, factors_t{f}), 1);
          if (! isempty (idx))
            factors_t{f} = dummy_bases{idx};
          endif
        endfor
        nm = strjoin (factors_t, ':');
        if (! any (strcmp (disp_terms, nm)))
          disp_terms{end+1} = nm;
        endif
      endfor

      if (has_intercept)
        lp_str = ifelse (isempty (disp_terms), '1', ['1 + ', strjoin(disp_terms, ' + ')]);
      else
        lp_str = strjoin (disp_terms, ' + ');
      endif

      FormulaS.ResponseName    = resp_name;
      FormulaS.LinearPredictor = lp_str;
      FormulaS.PredictorNames  = pred_names_raw;
      FormulaS.TermNames       = coef_names;
      FormulaS.HasIntercept    = has_intercept;
      FormulaS.Terms           = terms;
      FormulaS.InModel         = true (1, n_coef);
      FormulaS.NTerms          = n_coef;
      FormulaS.NPredictors     = p_raw;
      FormulaS.NVars           = n_vars;

      ObsInfo = table (w_full, excluded_mask, missing_mask, subset_mask, ...
        'VariableNames', {'Weights', 'Excluded', 'Missing', 'Subset'});

      if (! istable (data))
        VarsTable = array2table ([X_raw, y_full], 'VariableNames', var_names_all);
      else
        VarsTable = tbl;
      endif

      nv_total   = numel (var_names_all);
      vi_class   = cell  (nv_total, 1);
      vi_range   = cell  (nv_total, 1);
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
          fv = col_d(isfinite (col_d));
          vi_range{j} = ifelse (isempty (fv), [NaN, NaN], [min(fv), max(fv)]);
        else
          col_d = tbl.(vname);
          vi_class{j} = class (col_d);
          if (isnumeric (col_d))
            fv = col_d(isfinite (col_d));
            vi_range{j} = ifelse (isempty (fv), [NaN, NaN], [min(fv), max(fv)]);
          elseif (iscell (col_d))
            vi_range{j} = unique (col_d);
          else
            vi_range{j} = {};
          endif
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
      this.ModelFitVsNullModel      = struct ('Fstat',     Fstat, ...
                                             'Pvalue',    Fpval, ...
                                             'NullModel', 'constant');
      this.MSE                      = MSE;
      this.Residuals                = ResidTable;
      this.RMSE                     = RMSE;
      this.Rsquared                 = struct ('Ordinary', R2_ord, 'Adjusted', R2_adj);
      this.SSE                      = SSE;
      this.SSR                      = SSR;
      this.SST                      = SST;
      this.Robust                   = RobustS;
      this.Steps                    = [];
      this.Formula                  = FormulaS;
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
    ## @multitable @columnfractions 0.2 0.02 0.78
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"Alpha"} @tab @tab Significance level for the confidence
    ## interval, specified as a scalar in @math{[0,1]}.  The interval has
    ## coverage @math{100(1-\alpha)\%}.  Default is @code{0.05}, giving a 95%
    ## interval.
    ##
    ## @item @qcode{"Prediction"} @tab @tab Type of interval to compute.
    ## @code{"curve"} (default) gives a confidence interval on the mean response
    ## @math{f(x)}.  @code{"observation"} gives a wider prediction interval for
    ## a single future observation @math{y = f(x) + \varepsilon}, which accounts
    ## for both estimation uncertainty and irreducible noise; it adds
    ## @code{@var{mdl}.MSE} to the variance before computing the half-width.
    ##
    ## @item @qcode{"Simultaneous"} @tab @tab Logical flag controlling whether
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
            error ('predict: Alpha must be a scalar in [0,1].');
          endif
          i += 2;
        elseif (strcmpi (varargin{i}, 'Prediction'))
          pred_str = lower (char (varargin{i+1}));
          if (! any (strcmp (pred_str, {'curve', 'observation'})))
            error ('predict: Prediction must be ''curve'' or ''observation''.');
          endif
          pred_obs = strcmp (pred_str, 'observation');
          i += 2;
        elseif (strcmpi (varargin{i}, 'Simultaneous'))
          simultan = logical (varargin{i+1});
          i += 2;
        else
          error ('predict: unknown option ''%s''.', varargin{i});
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
            error ('predict: Xnew table is missing predictor ''%s''.', pred_names{j});
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
          error ('predict: Xnew must have %d columns.', p_raw);
        endif
        n_new = rows (X_raw);
      endif

      nan_rows     = any (isnan (X_raw), 2);
      X_enc_new    = LinearModel.lm_predict (X_raw, pred_names, mdl.CatLevelInfo, mdl.EncPredictorNames);
      X_design_new = LinearModel.lm_build_design (mdl.TermsMatrix, X_enc_new);

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
    ## @deftypefn {LinearModel} {@var{ysim} =} random (@var{mdl}, @var{Xnew})
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
        error ('random: Not enough input arguments.');
      endif
      if (nargin > 2)
        error ('random: Too many input arguments.');
      endif
      if (isempty (Xnew))
        error ('random: Xnew must have %d columns.', mdl.NumPredictors);
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
        error ('feval: Not enough input arguments.');
      endif

      if (n_extra == 1)

        Xnew = varargin{1};

        if (istable (Xnew))
          for j = 1:p_raw
            if (! ismember (mdl.PredictorNames{j}, Xnew.Properties.VariableNames))
              error (['feval: X does not contain one or more predictor ' ...
                      'variables needed for this model.']);
            endif
          endfor
        else
          if (columns (double (Xnew)) != p_raw)
            error ('feval: Predictor data matrix must have %d columns.', p_raw);
          endif
        endif

        ypred = predict (mdl, Xnew);

      elseif (n_extra == p_raw)

        ref_size = [];
        for i = 1:n_extra
          if (! isscalar (varargin{i}))
            s_i = size (varargin{i});
            if (isempty (ref_size))
              ref_size = s_i;
            elseif (! isequal (s_i, ref_size))
              error ('feval: All input arguments must be the same size.');
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

        error (['feval: Incorrect number of input arguments. You must provide ' ...
                'either %d separate predictor variable arguments, or one ' ...
                'predictor matrix with %d columns.'], p_raw, p_raw);

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
        error ('coefCI: Too many input arguments.');
      endif
      if (nargin < 2)
        alpha = 0.05;
      endif
      if (! isscalar (alpha))
        error (['coefCI: Invalid argument at position 2. ' ...
                'Value must be a scalar.']);
      endif
      if (! (alpha >= 0))
        error (['coefCI: Invalid argument at position 2. ' ...
                'Value must be greater than or equal to 0.']);
      endif
      if (alpha > 1)
        error (['coefCI: Invalid argument at position 2. ' ...
                'Value must be less than or equal to 1.']);
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
        error ('coefTest: Too many output arguments.');
      endif
      if (numel (varargin) > 2)
        error ('coefTest: Too many input arguments.');
      endif

      k = mdl.NumCoefficients;

      if (numel (varargin) >= 1 && ! isempty (varargin{1}))

        H = varargin{1};
        if (! isnumeric (H))
          error ('coefTest: H must be a %d-by-%d numeric matrix.', size (H, 1), k);
        endif
        if (size (H, 2) != k)
          error ('coefTest: H must be a %d-by-%d numeric matrix.', size (H, 1), k);
        endif
        if (any (any (isnan (H))))
          error (['coefTest: H is not full rank and hypotheses ' ...
                  'are not consistent.']);
        endif
        r = size (H, 1);

        if (numel (varargin) == 2)
          C = varargin{2};
          if (! isnumeric (C))
            error ('coefTest: C must be a numeric vector.');
          endif
          C = C(:);
          if (numel (C) != r)
            error ('coefTest: H must be a %d-by-%d numeric matrix.', numel (C), k);
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
        error ('dwtest: Too many output arguments.');
      endif
      if (numel (varargin) > 2)
        error ('dwtest: Too many input arguments.');
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
        error ('dwtest: The METHOD argument must be ''approximate'' or ''exact''.');
      endif
      method = lower (method);
      tail   = lower (tail);

      subset = logical (mdl.ObservationInfo.Subset);
      r      = mdl.Residuals.Raw(subset);
      r      = r(:);
      n      = numel (r);
      DW     = sum (diff (r).^2) / sum (r.^2);

      e      = ones (n-1, 1);
      A      = spdiags ([-e, e], [0, 1], n-1, n);
      M      = full (A' * A);
      [Q, R] = qr (mdl.DesignMatrix);
      dr     = abs (diag (R));
      tol    = max (size (mdl.DesignMatrix)) * eps (max (dr));
      rnk    = sum (dr > tol);
      Q2     = Q(:, rnk+1:size (Q, 2));
      lam    = real (eig (Q2' * M * Q2));
      k      = n - rnk;

      if (strcmp (method, 'exact'))
        a_row = (lam(:) - DW)';
        f     = @(u) sin (0.5 * sum (atan (u(:) * a_row), 2)) ./ ...
                     (u(:) .* prod ((1 + (u(:) * a_row).^2).^0.25, 2));
        p_r   = 0.5 - integral (f, 0, Inf, 'AbsTol', 1e-10, 'RelTol', 1e-8) / pi;
        p_r   = min (max (p_r, 0), 1);
      else
        mu_dw  = sum (lam) / k;
        var_dw = 2 * (sum (lam.^2) - sum (lam)^2 / k) / (k * (k + 2));
        p_r    = normcdf ((DW - mu_dw) / sqrt (var_dw));
        p_r    = min (max (p_r, 0), 1);
      endif

      if (strcmp (tail, 'right'))
        p = p_r;
      elseif (strcmp (tail, 'left'))
        p = 1 - p_r;
      else
        p = 2 * min (p_r, 1 - p_r);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {LinearModel} {@var{NewMdl} =} addTerms (@var{mdl}, @var{terms})
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
        error ('addTerms: Not enough input arguments.');
      endif
      if (nargin > 2)
        error ('addTerms: Too many input arguments.');
      endif

      nv   = mdl.NumVariables;
      pred = mdl.PredictorNames;

      if (isnumeric (terms) || islogical (terms))

        T = double (terms);
        if (isempty (T))
          error ('addTerms: Terms matrix must have %d columns.', nv);
        endif
        if (columns (T) == nv - 1)
          T = [T, zeros(rows (T), 1)];
        endif
        if (columns (T) != nv)
          error ('addTerms: Terms matrix must have %d columns.', nv);
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
                  error ('addTerms: Unrecognized variable: ''%s''.', vname);
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
                    error ('addTerms: Unrecognized variable: ''%s''.', vname);
                  endif
                  row(idx(1)) = row(idx(1)) + exp;
                endfor
                T = [T; row];
              endif
            endif
          endfor

        endfor

      else
        error (['addTerms: Model update specification must be a model formula', ...
                ' character vector or string scalar, or a terms matrix']);
      endif

      existing = mdl.TermsMatrix;
      is_new   = false (rows (T), 1);
      for i = 1:rows (T)
        is_new(i) = ! any (all (existing == T(i,:), 2));
      endfor
      new_rows = T(is_new, :);

      if (isempty (new_rows))
        warning ('addTerms: There are no new terms among the terms you specified.');
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
    ## @deftypefn {LinearModel} {@var{NewMdl} =} removeTerms (@var{mdl}, @var{terms})
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
        error ('removeTerms: Not enough input arguments.');
      endif
      if (nargin > 2)
        error ('removeTerms: Too many input arguments.');
      endif

      nv   = mdl.NumVariables;
      pred = mdl.PredictorNames;

      if (isnumeric (terms) || islogical (terms))

        T = double (terms);
        if (isempty (T))
          error ('removeTerms: Terms matrix must have %d columns.', nv);
        endif
        if (columns (T) == nv - 1)
          T = [T, zeros(rows (T), 1)];
        endif
        if (columns (T) != nv)
          error ('removeTerms: Terms matrix must have %d columns.', nv);
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
                  error ('removeTerms: Unrecognized variable: ''%s''.', vname);
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
                    error ('removeTerms: Unrecognized variable: ''%s''.', vname);
                  endif
                  row(idx(1)) = row(idx(1)) + exp;
                endfor
                T = [T; row];
              endif
            endif
          endfor

        endfor

      else
        error (['removeTerms: Model update specification must be a model formula', ...
                ' character vector or string scalar, or a terms matrix']);
      endif

      nc = columns (mdl.TermsMatrix);
      if (nc != nv)
        cat_info    = mdl.CatLevelInfo;
        n_pred      = nv - 1;
        orig_to_enc = cell (n_pred, 1);
        ecol        = 1;
        for j = 1:n_pred
          ci = [];
          if (! isempty (cat_info) && isfield (cat_info, 'names') ...
              && ! isempty (cat_info.names))
            ci = find (strcmp (cat_info.names, pred{j}));
          endif
          if (isempty (ci))
            orig_to_enc{j} = ecol;
            ecol            = ecol + 1;
          else
            n_lev           = numel (cat_info.levels{ci});
            orig_to_enc{j}  = ecol:(ecol + n_lev - 2);
            ecol            = ecol + n_lev - 1;
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
        warning ('removeTerms: No specified terms appear in the model.');
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
    ## @item @qcode{"histogram"} (default)
    ## Probability density histogram.  Only active observations are used.
    ## Returns one @code{patch} handle.  Accepts @code{FaceColor},
    ## @code{EdgeColor}, @code{FaceAlpha}, and @code{LineWidth} Name-Value
    ## arguments.
    ##
    ## @item @qcode{"fitted"}
    ## Residuals on the y-axis against fitted values on the x-axis.  A dotted
    ## horizontal reference line marks @math{y = 0}.  Returns two line handles:
    ## @code{h(1)} is the data scatter and @code{h(2)} is the reference line.
    ##
    ## @item @qcode{"caseorder"}
    ## Residuals on the y-axis against observation row number on the x-axis,
    ## covering all rows from 1 to @code{n_total}.  A dotted horizontal
    ## reference line marks @math{y = 0}.  Returns two line handles:
    ## @code{h(1)} is the data and @code{h(2)} is the reference line.
    ##
    ## @item @qcode{"lagged"}
    ## Each residual @math{r(t)} on the y-axis against the preceding residual
    ## @math{r(t-1)} on the x-axis.  Two dotted reference lines mark
    ## @math{y = 0} and @math{x = 0}.  Returns three line handles: @code{h(1)}
    ## is the scatter, @code{h(2)} is the horizontal reference, and @code{h(3)}
    ## is the vertical reference.
    ##
    ## @item @qcode{"probability"}
    ## Normal probability plot of the sorted active residuals produced by
    ## @code{normplot}.  Returns two handles: @code{h(1)} is the data line and
    ## @code{h(2)} is the fitted reference line produced by @code{normplot}.
    ## Name-Value arguments are not applied for this plot type.
    ##
    ## @item @qcode{"observed"}
    ## Observed response values on the y-axis against fitted values on the
    ## x-axis.  A dotted @math{y = x} reference line is drawn through the
    ## origin.  Vertical segments connect each observed point down to the
    ## reference line.
    ## Returns three handles: @code{h(1)} is the scatter, @code{h(2)} is the
    ## @math{y = x} reference, and @code{h(3)} is the vertical segment line
    ## (stored as a single @code{NaN}-separated line object).
    ##
    ## @item @qcode{"symmetry"}
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
    ## @multitable @columnfractions 0.28 0.02 0.70
    ## @headitem Name @tab @tab Description and default
    ##
    ## @item @qcode{"ResidualType"} @tab @tab
    ## Type of residual to plot.  One of @qcode{"raw"} (default),
    ## @qcode{"pearson"}, @qcode{"standardized"}, or @qcode{"studentized"}.
    ## Case-insensitive.  Selects the corresponding column of
    ## @code{mdl.Residuals}.
    ##
    ## @item @qcode{"Color"} @tab @tab
    ## (@emph{non-histogram}) Marker color.
    ## Default: @code{[0.1490 0.5490 0.8660]}.
    ##
    ## @item @qcode{"Marker"} @tab @tab
    ## (@emph{non-histogram}) Marker symbol.  Any symbol accepted by
    ## @code{plot} is valid.  Default: @qcode{"x"}.
    ##
    ## @item @qcode{"MarkerSize"} @tab @tab
    ## (@emph{non-histogram}) Marker size in points.  Default: @code{6}.
    ##
    ## @item @qcode{"MarkerEdgeColor"} @tab @tab
    ## (@emph{non-histogram}) Marker edge color.  Default: @qcode{"auto"}.
    ##
    ## @item @qcode{"MarkerFaceColor"} @tab @tab
    ## (@emph{non-histogram}) Marker fill color.  Default: @qcode{"none"}.
    ##
    ## @item @qcode{"LineWidth"} @tab @tab
    ## (@emph{non-histogram}) Width of the marker edge in points.
    ## Default: @code{0.5}.
    ##
    ## @item @qcode{"FaceColor"} @tab @tab
    ## (@emph{histogram only}) Fill color of the histogram bars.
    ## Default: @code{[0.1490 0.5490 0.8660]}.
    ##
    ## @item @qcode{"EdgeColor"} @tab @tab
    ## (@emph{histogram only}) Edge color of the histogram bars.
    ##
    ## @item @qcode{"FaceAlpha"} @tab @tab
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
          error ('plotResiduals: Bad residuals plot type.');
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
            error ('plotResiduals: ResidualType requires a value.');
          endif
          rt_val   = lower (char (args{i+1}));
          valid_rt = {'raw', 'pearson', 'standardized', 'studentized'};
          if (! any (strcmp (rt_val, valid_rt)))
            error (['plotResiduals: invalid ResidualType ''%s''. ', ...
                    'Valid values are: ''Raw'', ''Pearson'', ', ...
                    '''Standardized'', ''Studentized''.'], args{i+1});
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
          if     (frac <= 1); nice = 1;
          elseif (frac <= 2); nice = 2;
          elseif (frac <= 5); nice = 5;
          else;                nice = 10;
          endif
          bw = nice * mag;
          lo = floor (min (r_act) / bw) * bw;
          hi = ceil  (max (r_act) / bw) * bw;
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
          title  (ax, 'Histogram of residuals');

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
          title  (ax, 'Plot of residuals vs. fitted values');

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
          title  (ax, 'Case order plot of residuals');

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
          title  (ax, 'Plot of residuals vs. lagged residuals');

        case 'probability'
          r_act = r(! isnan (r));
          r_s   = sort (r_act);
          h     = normplot (ax, r_s);
          title  (ax, 'Normal probability plot of residuals');
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
          title  (ax, 'Plot of observed vs. fitted values');

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
          title  (ax, 'Symmetry plot of residuals around their median');

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
    ## @item @qcode{"leverage"} (default)
    ## Leverage of each observation (@code{mdl.Diagnostics.Leverage}).
    ## One dotted horizontal reference line at @math{2p/n}.
    ## Returns two handles: @code{h(1)} is the data scatter and @code{h(2)}
    ## is the reference line.
    ##
    ## @item @qcode{"cookd"}
    ## Cook's distance for each observation
    ## (@code{mdl.Diagnostics.CooksDistance}).  One dotted reference line at
    ## @math{3 \times \mathrm{mean(CooksDistance)}}, where the mean ignores
    ## @code{NaN} values.  Returns two handles: @code{h(1)} data, @code{h(2)}
    ## reference.
    ##
    ## @item @qcode{"covratio"}
    ## Delete-1 ratio of the determinant of the coefficient covariance matrix
    ## (@code{mdl.Diagnostics.CovRatio}).  Two dotted reference lines at
    ## @math{1 - 3p/n} (lower bound) and @math{1 + 3p/n} (upper bound).
    ## Both bounds are stored as a single @code{NaN}-separated line object.
    ## Returns two handles: @code{h(1)} data, @code{h(2)} combined reference.
    ##
    ## @item @qcode{"dfbetas"}
    ## Delete-1 scaled change in each coefficient estimate
    ## (@code{mdl.Diagnostics.Dfbetas}, one column per coefficient).
    ## One line object is drawn per coefficient.  Two dotted reference lines
    ## at @math{\pm 3/\sqrt{n}} are stored as a single @code{NaN}-separated
    ## line object.  Returns @math{p+1} handles: @code{h(1)} through
    ## @code{h(p)} are the per-coefficient data lines and @code{h(p+1)} is
    ## the combined reference.  Name-Value arguments are applied to all
    ## @math{p} data handles.
    ##
    ## @item @qcode{"dffits"}
    ## Delete-1 scaled change in the fitted value
    ## (@code{mdl.Diagnostics.Dffits}).  Two dotted reference lines at
    ## @math{\pm 2\sqrt{p/n}} stored as a single @code{NaN}-separated line.
    ## Returns two handles: @code{h(1)} data, @code{h(2)} combined reference.
    ##
    ## @item @qcode{"s2_i"}
    ## Delete-1 variance estimate (@code{mdl.Diagnostics.S2_i}).  One dotted
    ## reference line at @code{mdl.MSE}.  Returns two handles: @code{h(1)}
    ## data, @code{h(2)} reference.
    ##
    ## @item @qcode{"contour"}
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
    ## @multitable @columnfractions 0.28 0.02 0.70
    ## @headitem Name @tab @tab Description and default
    ##
    ## @item @qcode{"Color"} @tab @tab
    ## Marker color for data points.  For @code{"dfbetas"} this color is
    ## applied to all @math{p} coefficient line objects.
    ## Default: @code{[0.1490 0.5490 0.8660]}.
    ##
    ## @item @qcode{"Marker"} @tab @tab
    ## Marker symbol.  Any symbol accepted by @code{plot} is valid.
    ## Default: @qcode{"x"}.
    ##
    ## @item @qcode{"MarkerSize"} @tab @tab
    ## Marker size in points.  Default: @code{6}.
    ##
    ## @item @qcode{"MarkerEdgeColor"} @tab @tab
    ## Marker edge color.  Default: @qcode{"auto"}.
    ##
    ## @item @qcode{"MarkerFaceColor"} @tab @tab
    ## Marker fill color.  Default: @qcode{"none"}.
    ##
    ## @item @qcode{"LineWidth"} @tab @tab
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
          error ('plotDiagnostics: Bad diagnostics plot type.');
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
          title  (ax, 'Case order plot of leverage');

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
          title  (ax, 'Case order plot of Cook''s distance');

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
          title  (ax, 'Case order plot of covariance ratio');

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
          title  (ax, 'Case order plot of scaled change in coefficients');

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
          title  (ax, 'Case order plot of scaled change in fit');

        case 's2_i'
          s2 = diag_t.S2_i;
          hold (ax, 'on');
          h(1) = lm_plot_data (ax, 1:n, s2, props);
          h(2) = line ([0, n], [mse, mse], ...
                       'LineStyle', ':', 'Color', REF_COLOR, 'Parent', ax);
          hold (ax, 'off');
          xlabel (ax, 'Row number');
          ylabel (ax, 'Leave-one-out variance');
          title  (ax, 'Case order plot of leave-one-out variance');

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
          title  (ax, 'Cook''s distance factorization');

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
    ## @qcode{"varname: min to max"}, showing the predictor name and the
    ## minimum and maximum observed values used to compute the effect.
    ##
    ## @end deftypefn
    function h = plotEffects (this, varargin)
      [ax, mdl, args] = lm_plot_axes (this, varargin);

      if (! isempty (args))
        error ('plotEffects: Wrong number of arguments.');
      endif

      p = mdl.NumPredictors;
      if (! any (any (mdl.TermsMatrix(:, 1:end-1) != 0)))
        error ('plotEffects: Model has no predictors.');
      endif

      if (isempty (ax))
        ax = gca ();
      endif

      DEF_COLOR = [0.1490, 0.5490, 0.8660];

      act    = mdl.ObservationInfo.Subset;
      pred   = mdl.PredictorNames;
      V      = mdl.CoefficientCovariance;
      beta   = mdl.Coefficients.Estimate;
      t_crit = tinv (0.975, mdl.DFE);
      n_act  = sum (act);
      cinfo  = mdl.CatLevelInfo;
      ename  = mdl.EncPredictorNames;

      X_act    = zeros (n_act, p);
      is_cat   = false (1, p);
      cat_lvls = cell (1, p);

      for j = 1:p
        ci = [];
        if (! isempty (cinfo) && isfield (cinfo, 'names') ...
            && ! isempty (cinfo.names))
          ci = find (strcmp (cinfo.names, pred{j}));
        endif
        col = mdl.Variables{act, pred{j}};
        if (! isempty (ci))
          is_cat(j)   = true;
          levels_j    = cinfo.levels{ci};
          cat_lvls{j} = levels_j;
          if (iscell (col))
            col_str = col;
          elseif (isa (col, 'categorical'))
            col_str = cellstr (col);
          else
            col_str = cellstr (num2str (col(:)));
          endif
          codes = zeros (n_act, 1);
          for L = 1:numel (levels_j)
            codes (strcmp (col_str, char (levels_j{L}))) = L;
          endfor
          X_act(:,j) = codes;
        else
          X_act(:,j) = double (col(:));
        endif
      endfor

      effects = zeros (1, p);
      ci_lo   = zeros (1, p);
      ci_hi   = zeros (1, p);
      x_lo_v  = zeros (1, p);
      x_hi_v  = zeros (1, p);

      for j = 1:p
        if (is_cat(j))
          x_lo_v(j) = 1;
          x_hi_v(j) = numel (cat_lvls{j});
        else
          x_lo_v(j) = min (X_act(:,j));
          x_hi_v(j) = max (X_act(:,j));
        endif

        X_hi_rows = X_act;  X_hi_rows(:,j) = x_hi_v(j);
        X_lo_rows = X_act;  X_lo_rows(:,j) = x_lo_v(j);

        X_hi_enc = LinearModel.lm_predict (X_hi_rows, pred, cinfo, ename);
        X_lo_enc = LinearModel.lm_predict (X_lo_rows, pred, cinfo, ename);
        D_hi     = LinearModel.lm_build_design (mdl.TermsMatrix, X_hi_enc);
        D_lo     = LinearModel.lm_build_design (mdl.TermsMatrix, X_lo_enc);

        c_bar      = mean (D_hi - D_lo, 1);
        effects(j) = c_bar * beta;
        SE         = sqrt (max (0, c_bar * V * c_bar'));

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

      ytl = cell (p, 1);
      for j = 1:p
        if (is_cat(j))
          lo_str = char (cat_lvls{j}{x_lo_v(j)});
          hi_str = char (cat_lvls{j}{x_hi_v(j)});
        else
          lo_str = num2str (x_lo_v(j), '%g');
          hi_str = num2str (x_hi_v(j), '%g');
        endif
        ytl{j} = [pred{j}, ': ', lo_str, ' to ', hi_str];
      endfor

      set (ax, 'YTick', 1:p, 'YTickLabel', ytl, 'YDir', 'reverse');
      ylim  (ax, [0.5, p + 0.5]);
      xlabel (ax, 'Main Effect');
      ylabel (ax, '');
      title  (ax, 'Main Effects Plot');

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
    ## @multitable @columnfractions 0.28 0.02 0.70
    ## @headitem Name @tab @tab Description and default
    ##
    ## @item @qcode{"Color"} @tab @tab
    ## Marker color.  Default: @code{[0.1490 0.5490 0.8660]}.
    ##
    ## @item @qcode{"Marker"} @tab @tab
    ## Marker symbol.  Default: @qcode{"x"}.
    ##
    ## @item @qcode{"MarkerSize"} @tab @tab
    ## Marker size in points.  Default: @code{6}.
    ##
    ## @item @qcode{"MarkerEdgeColor"} @tab @tab
    ## Marker edge color.  Default: @qcode{"auto"}.
    ##
    ## @item @qcode{"MarkerFaceColor"} @tab @tab
    ## Marker fill color.  Default: @qcode{"none"}.
    ##
    ## @item @qcode{"LineWidth"} @tab @tab
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
        error ('plotAdjustedResponse: Not enough input arguments.');
      endif

      var  = args{1};
      args = args(2:end);

      vnames = mdl.VariableNames;
      pred   = mdl.PredictorNames;

      if (ischar (var) || isstring (var))
        vname = char (var);
        if (isempty (find (strcmp (vnames, vname))))
          error ('plotAdjustedResponse: ''%s'' is not a variable for this fit.', vname);
        endif
      elseif (isnumeric (var) && isscalar (var))
        if (var != fix (var) || var < 1)
          error (['plotAdjustedResponse: Variable must be specified as a ', ...
                  'name or a positive integer.']);
        endif
        if (var > numel (vnames))
          error ('plotAdjustedResponse: This model only contains %d variables.', ...
                 numel (vnames));
        endif
        vname = vnames{var};
      else
        error (['plotAdjustedResponse: Variable must be specified as a ', ...
                'name or a positive integer.']);
      endif

      if (strcmp (vname, mdl.ResponseName))
        error ('plotAdjustedResponse: The variable ''%s'' is the response in this model.', ...
               vname);
      endif

      j     = find (strcmp (pred, vname));
      pname = pred{j};

      act   = mdl.ObservationInfo.Subset;
      cinfo = mdl.CatLevelInfo;
      n_act = sum (act);
      p     = numel (pred);

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
          if (iscell (col))
            col_str = col;
          elseif (isa (col, 'categorical'))
            col_str = cellstr (col);
          else
            col_str = cellstr (num2str (col(:)));
          endif
          codes = zeros (n_act, 1);
          for L = 1:numel (levels_k)
            codes (strcmp (col_str, char (levels_k{L}))) = L;
          endfor
          X_act(:,k) = codes;
        else
          X_act(:,k) = double (col(:));
        endif
      endfor

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
          X_enc = LinearModel.lm_predict (X_rows, pred, mdl.CatLevelInfo, mdl.EncPredictorNames);
          D     = LinearModel.lm_build_design (mdl.TermsMatrix, X_enc);
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
          X_enc = LinearModel.lm_predict (X_rows, pred, mdl.CatLevelInfo, mdl.EncPredictorNames);
          D     = LinearModel.lm_build_design (mdl.TermsMatrix, X_enc);
          g_active(i) = mean (D * beta);
        endfor
        ydata(act) = g_active + resid(act);

        fit_x = linspace (min (x_active), max (x_active), 100)';
        fit_y = zeros (100, 1);
        for k = 1:100
          X_rows      = X_act;
          X_rows(:,j) = fit_x(k);
          X_enc = LinearModel.lm_predict (X_rows, pred, mdl.CatLevelInfo, mdl.EncPredictorNames);
          D     = LinearModel.lm_build_design (mdl.TermsMatrix, X_enc);
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
      title  (ax, 'Adjusted Response Plot');
      
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
    ## @multitable @columnfractions 0.28 0.02 0.70
    ## @headitem Name @tab @tab Description and default
    ##
    ## @item @qcode{"Color"} @tab @tab
    ## Marker color.  Default: @code{[0.1490 0.5490 0.8660]}.
    ##
    ## @item @qcode{"Marker"} @tab @tab
    ## Marker symbol.  Default: @qcode{"x"}.
    ##
    ## @item @qcode{"MarkerSize"} @tab @tab
    ## Marker size in points.  Default: @code{6}.
    ##
    ## @item @qcode{"MarkerEdgeColor"} @tab @tab
    ## Marker edge color.  Default: @qcode{"auto"}.
    ##
    ## @item @qcode{"MarkerFaceColor"} @tab @tab
    ## Marker fill color.  Default: @qcode{"none"}.
    ##
    ## @item @qcode{"LineWidth"} @tab @tab
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
        label = 'Whole Model';
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
              error ('plotAdded: Bad coefficient name.');
            endif
            levels_ci = cinfo.levels{ci};
            J = zeros (1, numel (levels_ci) - 1);
            for L = 2:numel (levels_ci)
              cn = sprintf ('%s_%s', cname, char (levels_ci{L}));
              J(L-1) = find (strcmp (cnames, cn));
            endfor
          endif
        else
          J = coefarg(:)';
          if (any (J != fix (J)) || any (J < 1) || any (J > ncoef))
            error ('plotAdded: Bad coefficient number.');
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
                cn = sprintf ('%s_%s', cinfo.names{ci}, char (levels_ci{L}));
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
        error ('plotAdded: Bad coefficient number.');
      endif

      act   = mdl.ObservationInfo.Subset;
      cinfo = mdl.CatLevelInfo;
      n_act = sum (act);
      p     = numel (pred);

      X_act = zeros (n_act, p);
      for k = 1:p
        ci = [];
        if (! isempty (cinfo) && isfield (cinfo, 'names') ...
            && ! isempty (cinfo.names))
          ci = find (strcmp (cinfo.names, pred{k}));
        endif
        col = mdl.Variables{act, pred{k}};
        if (! isempty (ci))
          levels_k = cinfo.levels{ci};
          if (iscell (col))
            col_str = col;
          elseif (isa (col, 'categorical'))
            col_str = cellstr (col);
          else
            col_str = cellstr (num2str (col(:)));
          endif
          codes = zeros (n_act, 1);
          for L = 1:numel (levels_k)
            codes (strcmp (col_str, char (levels_k{L}))) = L;
          endfor
          X_act(:,k) = codes;
        else
          X_act(:,k) = double (col(:));
        endif
      endfor

      D = LinearModel.lm_predict (X_act, pred, cinfo, mdl.EncPredictorNames);
      D = LinearModel.lm_build_design (mdl.TermsMatrix, D);

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
                   'DisplayName', sprintf ('Fit: y = %g*x', slope));
      h(3) = line (bound_x, bound_y, 'Color', FIT_COLOR, 'LineStyle', ':', ...
                   'Marker', 'none', 'Parent', ax, ...
                   'DisplayName', '95% conf. bounds');
      hold (ax, 'off');

      xlabel (ax, ['Adjusted ', label]);
      ylabel (ax, ['Adjusted ', mdl.ResponseName]);
      title  (ax, ['Added Variable Plot for ', label]);
      
      hleg = legend (ax, 'show');
      set (hleg, 'Location', lm_legend_corner (xdata, ydata));

      if (nargout == 0)
        clear h;
      endif

    endfunction

  endmethods

  methods (Access = private, Static)

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

      n_eff = sum (w > 0);
      SSE   = sum (w .* Raw.^2);
      wmean = sum (w .* y) / max (sum (w), eps);
      SST   = sum (w .* (y - wmean).^2);
      SSR   = SST - SSE;

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
    endfunction

    function X_design = lm_build_design (terms, X_enc)
      n_obs    = rows (X_enc);
      n_coef   = rows (terms);
      p_enc    = columns (X_enc);
      X_design = zeros (n_obs, n_coef);
      for t = 1:n_coef
        term_row = terms(t, 1:p_enc);
        col_t    = ones (n_obs, 1);
        for j = find (term_row != 0)
          col_t = col_t .* (X_enc(:, j) .^ term_row(j));
        endfor
        X_design(:, t) = col_t;
      endfor
    endfunction

    function crit = lm_criteria (fit, n_obs, has_intercept)
      p   = fit.rank_X;
      SSE = fit.SSE;
      SSR = fit.SSR;
      SST = fit.SST;
      DFE = fit.DFE;
      MSE = fit.MSE;

      LogLikelihood = -(n_obs / 2) * (1 + log (2 * pi * SSE / n_obs));

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

    function X_enc = lm_predict (X_raw, pred_names, cat_info, enc_names)
      if (nargin < 4)
        enc_names = pred_names;
      endif
      n     = rows (X_raw);
      X_enc = zeros (n, numel (enc_names));
      for c = 1:numel (enc_names)
        name = enc_names{c};

        j = find (strcmp (pred_names, name), 1);
        if (! isempty (j))
          X_enc(:, c) = X_raw(:, j);
          continue;
        endif

        tok = regexp (name, '^(.+)\^(\d+)$', 'tokens');
        if (! isempty (tok))
          j = find (strcmp (pred_names, tok{1}{1}), 1);
          k = str2double (tok{1}{2});
          X_enc(:, c) = X_raw(:, j) .^ k;
          continue;
        endif

        found = false;
        for j = 1:numel (pred_names)
          ci = [];
          if (! isempty (cat_info.names))
            ci = find (strcmp (cat_info.names, pred_names{j}));
          endif
          if (isempty (ci))
            continue;
          endif
          levels_j = cat_info.levels{ci};
          for L = 2:numel (levels_j)
            if (strcmp (name, sprintf ('%s_%s', pred_names{j}, char (levels_j{L}))))
              X_enc(:, c) = double (X_raw(:, j) == L);
              found = true;
              break;
            endif
          endfor
          if (found)
            break;
          endif
        endfor
      endfor
    endfunction

  endmethods

endclassdef

function opts = lm_parse_nv (nv_args)

  opt_names = {'Intercept', 'Weights', 'Exclude', 'RobustOpts', ...
               'VarNames', 'CategoricalVars', 'ResponseVar', 'PredictorVars'};
  def_vals  = {true, [], [], [], {}, [], '', {}};
  [intercept, weights, exclude, robustopts, varnames, catvars, ...
   respvar, predvars, rem_args] = parsePairedArguments (opt_names, def_vals, nv_args);

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
  else
    opts.PredictorVars = cellstr (predvars);
  endif

endfunction

## Parse modelspec into terms matrix.
function [terms, has_intercept, coef_names] = lm_parse_modelspec ( ...
    modelspec, pred_names, n_preds, intercept_nv)

  p = n_preds;

  if (isempty (modelspec) || (ischar (modelspec) && strcmpi (modelspec, 'linear')))
    terms = [zeros(1, p+1); [eye(p), zeros(p, 1)]];

  elseif (ischar (modelspec) && strcmpi (modelspec, 'constant'))
    terms = zeros (1, p+1);

  elseif (ischar (modelspec) && strcmpi (modelspec, 'interactions'))
    linear_part = [zeros(1, p+1); [eye(p), zeros(p, 1)]];
    inter_part  = zeros (0, p+1);
    for i = 1:p
      for j = i+1:p
        row = zeros (1, p+1); row(i) = 1; row(j) = 1;
        inter_part = [inter_part; row];
      endfor
    endfor
    terms = [linear_part; inter_part];

  elseif (ischar (modelspec) && strcmpi (modelspec, 'purequadratic'))
    linear_part = [zeros(1, p+1); [eye(p), zeros(p, 1)]];
    quad_part   = zeros (p, p+1);
    for j = 1:p
      quad_part(j, j) = 2;
    endfor
    terms = [linear_part; quad_part];

  elseif (ischar (modelspec) && strcmpi (modelspec, 'quadratic'))
    linear_part = [zeros(1, p+1); [eye(p), zeros(p, 1)]];
    quad_part   = zeros (p, p+1);
    for j = 1:p
      quad_part(j, j) = 2;
    endfor
    inter_part = zeros (0, p+1);
    for i = 1:p
      for j = i+1:p
        row = zeros (1, p+1); row(i) = 1; row(j) = 1;
        inter_part = [inter_part; row];
      endfor
    endfor
    terms = [linear_part; inter_part; quad_part];

  elseif (ischar (modelspec) && strcmpi (modelspec, 'full'))
    terms = zeros (1, p+1);
    for k = 1:p
      idx_mat = nchoosek (1:p, k);
      for j = 1:rows (idx_mat)
        row               = zeros (1, p+1);
        row(idx_mat(j,:)) = 1;
        terms             = [terms; row];
      endfor
    endfor

  elseif (isnumeric (modelspec))
    terms = double (modelspec);
    if (size (terms, 2) == p)
      terms = [terms, zeros(rows (terms), 1)];
    elseif (size (terms, 2) == p + 1)
      if (! all (terms(:, end) == 0))
        error ('LinearModel: Last column of terms matrix must be all zeros.');
      endif
    else
      error ('LinearModel: Terms matrix must have %d or %d columns.', p, p+1);
    endif

  else
    error ('fitlm: Unknown model specification.');
  endif

  if (! intercept_nv)
    int_rows = all (terms(:, 1:end-1) == 0, 2);
    terms    = terms(! int_rows, :);
  endif

  has_intercept = any (all (terms(:, 1:end-1) == 0, 2));

  n_terms    = rows (terms);
  coef_names = cell (1, n_terms);
  for t = 1:n_terms
    term_row = terms(t, 1:end-1);
    if (all (term_row == 0))
      coef_names{t} = '(Intercept)';
    else
      parts_t = {};
      for j = 1:numel (term_row)
        if (term_row(j) != 0)
          if (term_row(j) == 1)
            parts_t{end+1} = pred_names{j};
          else
            parts_t{end+1} = sprintf ('%s^%d', pred_names{j}, term_row(j));
          endif
        endif
      endfor
      coef_names{t} = strjoin (parts_t, ':');
    endif
  endfor
endfunction

## expand categorical columns into L-1 dummy variables
function [X_enc, enc_names, cat_info] = lm_encode_categorical ( ...
    X_num, cat_cols, pred_names, cat_levels)

  X_enc     = zeros (rows (X_num), 0);
  enc_names = {};
  cat_info.names  = {};
  cat_info.levels = {};

  for j = 1:numel (pred_names)
    if (! cat_cols(j))
      X_enc     = [X_enc, X_num(:, j)];
      enc_names = [enc_names, pred_names{j}];
    else
      levels_j = cat_levels{j};
      if (isempty (levels_j))
        uvals    = sort (unique (X_num(isfinite (X_num(:,j)), j)));
        levels_j = cellstr (num2str (uvals(:)));
      endif
      n_lev = numel (levels_j);
      for L = 2:n_lev
        dummy     = double (X_num(:, j) == L);
        X_enc     = [X_enc, dummy];
        enc_names = [enc_names, [pred_names{j}, '_', char(levels_j{L})]];
      endfor
      cat_info.names{end+1}  = pred_names{j};
      cat_info.levels{end+1} = levels_j;
    endif
  endfor
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
    error ('lm_plot_props: unrecognized property ''%s''.', rem_args{1});
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

function fit = lm_robust_fit (X, y, w, wgtfun, tune)

  n  = rows (X);
  p  = columns (X);
  w  = w(:);
  sw = sqrt (w);

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

    wts   = wgtfun (radj / (max (s, tiny_s) * tune));
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
    phi = u .* wgtfun (u);

    delta = 0.0001;
    u1   = u - delta;
    phi0 = u1 .* wgtfun (u1);
    u1   = u + delta;
    phi1 = u1 .* wgtfun (u1);
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
%! assert (mdl.NumObservations,          20);
%! assert (mdl.NumCoefficients,           3);
%! assert (mdl.NumVariables,              3);
%! assert (mdl.NumPredictors,             2);
%! assert (mdl.NumEstimatedCoefficients,  3);
%! assert (mdl.DFE,                      17);
%! assert (mdl.SSE,  0.386545331386823,   1e-9);
%! assert (mdl.SSR,  583.523874670959,    1e-6);
%! assert (mdl.SST,  583.910420002346,    1e-6);
%! assert (mdl.MSE,  0.0227379606698351,  1e-10);
%! assert (mdl.RMSE, 0.150791116017606,   1e-10);
%! assert (mdl.Rsquared.Ordinary, 0.999338005765704, 1e-10);
%! assert (mdl.Rsquared.Adjusted, 0.999260124091081, 1e-10);
%! assert (mdl.LogLikelihood, 11.0836133807695, 1e-6);
%! assert (mdl.ModelCriterion.AIC,  -16.1672267615389, 1e-6);
%! assert (mdl.ModelCriterion.AICc, -14.6672267615389, 1e-6);
%! assert (mdl.ModelCriterion.BIC,  -13.180029940877,  1e-6);
%! assert (mdl.ModelCriterion.CAIC, -10.180029940877,  1e-6);
%! assert (mdl.ModelFitVsNullModel.Fstat, 12831.4909842738, 1e-4);
%! assert (strcmp (mdl.ModelFitVsNullModel.NullModel, 'constant'));

%!test
%! ## constant-only model: SSR is exactly zero, SSE equals SST
%! mc = fitlm (X, y, 'constant');
%! assert (mc.SSE, 583.910420002346, 1e-6);
%! assert (mc.SSR, 0, 1e-12);
%! assert (mc.SSE, mc.SST, 1e-12);

%!test
%! ## coefficient estimates, SE, tStat, names, covariance, schema
%! assert (mdl.Coefficients.Estimate, [0.1161886778; 2.508451491; -0.9788353298], 1e-7);
%! assert (mdl.Coefficients.SE,       [0.112185831;  0.4920818186; 0.02276108523], 1e-8);
%! assert (mdl.Coefficients.tStat,    [1.035680502;  5.097630913; -43.00477415],   1e-6);
%! assert (all (mdl.Coefficients.pValue >= 0 & mdl.Coefficients.pValue <= 1));
%! assert (isequal (mdl.CoefficientNames, {'(Intercept)', 'x1', 'x2'}));
%! assert (isequal (mdl.CoefficientNames, mdl.Coefficients.Properties.RowNames(:)'));
%! assert (size (mdl.CoefficientCovariance), [3, 3]);
%! assert (diag (mdl.CoefficientCovariance), [0.0125857; 0.242145; 0.000518067], 1e-6);
%! assert (width (mdl.Coefficients), 4);
%! assert (isequal (mdl.Coefficients.Properties.VariableNames, ...
%!                  {'Estimate','SE','tStat','pValue'}));

%!test
%! ## fitted values, predict(), residual columns (obs 1-3), schema
%! assert (mdl.Fitted, y - mdl.Residuals.Raw, 1e-10);
%! yp = predict (mdl, X);
%! assert (size (yp), [20, 1]);
%! assert (yp(1), 0.192669485827491, 1e-10);
%! assert (yp(2), 0.171266760882256, 1e-10);
%! assert (mdl.Residuals.Raw(1:3),          [0.075624711134088; 0.110592724482880; -0.023756501342530], 1e-10);
%! assert (mdl.Residuals.Pearson(1:3),      [0.501519672586403; 0.733416711830473; -0.157545762442370], 1e-9);
%! assert (mdl.Residuals.Standardized(1:3), [0.632246516521578; 0.844226394951239; -0.172381368754725], 1e-8);
%! assert (mdl.Residuals.Studentized(1:3),  [0.620710275056923; 0.836747864205268; -0.167380843634378], 1e-6);
%! assert (width (mdl.Residuals), 4);
%! assert (isequal (mdl.Residuals.Properties.VariableNames, ...
%!                  {'Raw','Pearson','Studentized','Standardized'}));

%!test
%! ## diagnostics
%! H = mdl.Diagnostics.HatMatrix;
%! assert (size (H), [20, 20]);
%! assert (H, H', 1e-10);
%! assert (H * H, H, 1e-8);
%! assert (H(1,1), 0.370779220779221, 1e-10);
%! assert (H(1,2), 0.298051948051948, 1e-10);
%! assert (mdl.Diagnostics.Leverage(1:3),      [0.370779220779221; 0.245283663704716; 0.164718614718615], 1e-10);
%! assert (mdl.Diagnostics.CooksDistance(1:3), [0.078517048682575; 0.077211407930332; 0.001953301452841], 1e-8);
%! assert (mdl.Diagnostics.S2_i(1:3),          [0.023591009857798; 0.023146223303229; 0.024116854077430], 1e-8);
%! assert (mdl.Diagnostics.CovRatio(1:3),      [1.774933176573401; 1.397661919176034; 1.428481535363283], 1e-6);
%! assert (mdl.Diagnostics.Dffits(1:3),        [0.476480465355394; 0.477020506700835; -0.074329411030064], 1e-6);
%! assert (size (mdl.Diagnostics.Dfbetas), [20, 3]);
%! assert (width (mdl.Diagnostics), 7);
%! assert (isequal (mdl.Diagnostics.Properties.VariableNames, ...
%!                  {'Leverage','CooksDistance','Dffits','S2_i', ...
%!                   'CovRatio','Dfbetas','HatMatrix'}));

%!test
%! ## ObservationInfo, VariableInfo, names, Formula, Variables
%! assert (width (mdl.ObservationInfo), 4);
%! assert (height (mdl.ObservationInfo), 20);
%! assert (isequal (mdl.ObservationInfo.Properties.VariableNames, ...
%!                  {'Weights','Excluded','Missing','Subset'}));
%! assert (all (mdl.ObservationInfo.Weights == 1));
%! assert (all (mdl.ObservationInfo.Subset == ...
%!              (! mdl.ObservationInfo.Missing & ! mdl.ObservationInfo.Excluded)));
%! assert (width  (mdl.VariableInfo), 4);
%! assert (height (mdl.VariableInfo), 3);
%! assert (isequal (mdl.VariableInfo.Properties.VariableNames, ...
%!                  {'Class','Range','InModel','IsCategorical'}));
%! assert (mdl.VariableInfo.InModel(strcmp (mdl.VariableNames, 'y')), false);
%! assert (all (mdl.VariableInfo.InModel(! strcmp (mdl.VariableNames, 'y'))));
%! assert (mdl.ResponseName, 'y');
%! assert (isequal (mdl.PredictorNames, {'x1','x2'}));
%! assert (isequal (mdl.VariableNames, {'x1','x2','y'}));
%! assert (mdl.Formula.HasIntercept, true);
%! assert (mdl.Formula.LinearPredictor, '1 + x1 + x2');
%! assert (mdl.Formula.NTerms, 3);
%! assert (strcmp (mdl.Variables.Properties.VariableNames{end}, 'y'));

%!test
%! ## NaN in predictor drops the row from the fit
%! X2 = X;  X2(2,1) = NaN;
%! m = fitlm (X2, y);
%! assert (m.NumObservations, 19);
%! assert (m.ObservationInfo.Missing(2), true);
%! assert (m.ObservationInfo.Subset(2),  false);
%! assert (isnan (m.Fitted(2)));
%! assert (m.SSE, 0.370339572851658, 1e-9);
%! assert (m.SST, 547.616796178045,  1e-6);
%! assert (m.Coefficients.Estimate, ...
%!         [0.0641300185953764; 2.68263140079657; -0.985345792254554], 1e-7);
%! yp = predict (m, X2);
%! assert (isnan (yp(2)));
%! assert (! isnan (yp(1)));
%! assert (size (m.Diagnostics.HatMatrix), [20, 20]);
%! assert (m.Diagnostics.Leverage(1), 0.488485648300892, 1e-8);
%! assert (m.Diagnostics.CooksDistance(1), 0.38266162627176, 1e-6);

%!test
%! ## NaN in response drops the row but predict still works normally since X has no NaN
%! y3 = y;  y3(5) = NaN;
%! m = fitlm (X, y3);
%! assert (m.NumObservations, 19);
%! assert (m.ObservationInfo.Missing(5), true);
%! assert (isnan (m.Fitted(5)));
%! assert (m.SSE, 0.337042910721425, 1e-9);
%! assert (m.SST, 558.654961265991,  1e-6);
%! assert (m.Coefficients.Estimate, ...
%!         [0.145131680993155; 2.4865383021829; -0.979635081226208], 1e-7);
%! yp = predict (m, X);
%! assert (yp(5), -0.45777759499388, 1e-8);
%! assert (yp(1), 0.22047684204099, 1e-8);
%! assert (size (m.Diagnostics.HatMatrix), [20, 20]);
%! assert (m.Diagnostics.Leverage(1), 0.386399650026734, 1e-8);

%!test
%! ## multiple NaN rows drop all affected observations from the fit
%! X4 = X;  X4([2,8,14],2) = NaN;
%! m = fitlm (X4, y);
%! assert (sum (m.ObservationInfo.Missing), 3);
%! assert (m.NumObservations, 17);
%! assert (m.SSE, 0.261285495635633, 1e-9);
%! assert (m.SSR, 527.635694805749,  1e-6);
%! assert (m.SST, 527.896980301385,  1e-6);
%! assert (m.Coefficients.Estimate, ...
%!         [0.0986395043600395; 2.3735792982821; -0.97106191310122], 1e-7);
%! assert (size (m.Diagnostics.HatMatrix), [20, 20]);
%! assert (sum (m.Diagnostics.Leverage), 3, 1e-8);

%!test
%! ## exclude by index and exclude by logical vector give identical results
%! m  = fitlm (X, y, 'Exclude', [3, 7]);
%! excl = false (n, 1);  excl([3, 7]) = true;
%! m2 = fitlm (X, y, 'Exclude', excl);
%! assert (m.NumObservations, 18);
%! assert (sum (m.ObservationInfo.Excluded), 2);
%! assert (isnan (m.Fitted(3)) && isnan (m.Fitted(7)));
%! assert (m.Coefficients.Estimate, m2.Coefficients.Estimate, 1e-12);
%! assert (m.Coefficients.Estimate, ...
%!         [0.118938102486219; 2.43606890944554; -0.974833228191174], 1e-7);
%! ype = predict (m);
%! assert (size (ype), [20, 1]);
%! assert (! isnan (ype(3)) && ! isnan (ype(7)));
%! [~, yci] = predict (m);
%! assert (yci(1,1), -0.0283122762458446, 1e-10);
%! assert (yci(1,2),  0.412312049343719,  1e-10);
%! assert (size (m.Diagnostics.HatMatrix), [20, 20]);
%! assert (m.Diagnostics.Leverage(1), 0.437780279893411, 1e-8);
%! assert (m.Diagnostics.CooksDistance(1), 0.110112457355807, 1e-7);

%!test
%! ## NaN and exclude together remove both the missing and the excluded row
%! X6 = X;  X6(1,1) = NaN;
%! m = fitlm (X6, y, 'Exclude', [2]);
%! assert (m.NumObservations, 18);
%! assert (m.ObservationInfo.Missing(1),  true);
%! assert (m.ObservationInfo.Excluded(2), true);
%! assert (m.SSE, 0.342515396265007, 1e-9);
%! assert (m.Coefficients.Estimate, ...
%!         [-0.0735450184226009; 3.17679029176988; -1.0045827469016], 1e-7);

%!test
%! ## weighted least squares produces different SSE and stores the weights
%! w = abs (sin ((1:n)')) + 0.1;
%! m = fitlm (X, y, 'Weights', w);
%! assert (m.SSE, 0.363519720897775, 1e-10);
%! assert (m.ObservationInfo.Weights, w, 1e-15);
%! assert (m.SST, 4.419834786423099e+02, 1e-8);
%! [yp, yci] = predict (m, [0.5 0.25; 1.0 1.0]);
%! assert (yp(1),  1.106748776307639, 1e-10);
%! assert (yp(2),  1.593185531572655, 1e-10);
%! assert (yci(1,1), 0.763985050242272, 1e-10);
%! assert (yci(1,2), 1.449512502373006, 1e-10);
%! assert (m.Diagnostics.Leverage(1), 0.421642939812731, 1e-8);
%! assert (m.Diagnostics.Leverage(2), 0.301314342928707, 1e-8);
%! assert (m.Diagnostics.HatMatrix(1,1), 0.421642939812731, 1e-8);
%! assert (m.Diagnostics.CooksDistance(1), 0.0728569335883748, 1e-7);
%! assert (m.Diagnostics.CovRatio(1), 1.96611264276187, 1e-6);

%!test
%! ## uniform weights scale internals but leave point estimates unchanged
%! m = fitlm (X, y, 'Weights', 2 * ones (n, 1));
%! assert (m.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-10);

%!test
%! ## constant linear and default modelspecs behave as expected
%! m  = fitlm (X, y, 'constant');
%! assert (m.NumCoefficients, 1);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! m2 = fitlm (X, y, 'linear');
%! m3 = fitlm (X, y, []);
%! assert (m2.NumCoefficients, 3);
%! assert (m2.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-12);
%! assert (m3.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-12);

%!test
%! ## purequadratic modelspec produces the expected term count
%! m = fitlm (X, y, 'purequadratic');
%! assert (m.NumCoefficients, 5);

%!test
%! ## interactions modelspec term count and coefficients are verified
%! m = fitlm (X, y, 'interactions');
%! assert (m.NumCoefficients, 4);
%! assert (m.SSE, 0.383859187927621, 1e-9);
%! assert (m.Coefficients.Estimate, ...
%!         [0.157640728038039; 2.08542680311791; -0.929682701072813; -0.031208018255475], 1e-7);

%!test
%! ## quadratic modelspec is rank deficient for this design and drops one coefficient
%! m = fitlm (X, y, 'quadratic');
%! assert (m.NumCoefficients, 6);
%! assert (m.SSE, 0.315784637443501, 1e-9);
%! assert (m.Coefficients.Estimate, ...
%!         [0.447436249544699; -2.44859403731902; -0.0121968798776254; ...
%!          -1.36755100280532; 0; 0.0318176901083297], 1e-7);
%! drop = find (m.Coefficients.SE == 0);
%! assert (numel (drop), 1);
%! assert (isnan (m.Coefficients.tStat(drop)));

%!test
%! ## full modelspec with two predictors matches interactions exactly
%! m  = fitlm (X, y, 'full');
%! m2 = fitlm (X, y, 'interactions');
%! assert (m.NumCoefficients, 4);
%! assert (m.Coefficients.Estimate, m2.Coefficients.Estimate, 1e-10);
%! assert (m.Coefficients.Estimate, ...
%!         [0.157640728038039; 2.08542680311791; -0.929682701072813; -0.031208018255475], 1e-7);

%!test
%! ## full modelspec with three predictors includes the three way interaction term
%! X3 = [X, cos((1:n)' * pi / n)];
%! m  = fitlm (X3, y, 'full');
%! assert (m.NumCoefficients, 8);
%! assert (any (strcmp (m.CoefficientNames, 'x1:x2:x3')));
%! assert (m.SSE, 0.231331066631196, 1e-8);
%! idx3 = find (strcmp (m.CoefficientNames, 'x1:x2:x3'));
%! assert (m.Coefficients.Estimate(idx3), 0.514890561912964, 1e-6);

%!test
%! ## full modelspec without an intercept drops the intercept coefficient
%! m = fitlm (X, y, 'full', 'Intercept', false);
%! assert (m.NumCoefficients, 3);
%! assert (! any (strcmp (m.CoefficientNames, '(Intercept)')));
%! assert (m.Coefficients.Estimate, ...
%!         [3.232987312533958; -1.041484635851565; 0.0324190990982863], 1e-7);

%!test
%! ## a p column terms matrix produces a model with no intercept
%! m = fitlm (X, y, [1 0; 0 1]);
%! assert (m.NumCoefficients, 2);
%! assert (! any (strcmp (m.CoefficientNames, '(Intercept)')));
%! assert (m.Coefficients.Estimate, [2.96142161317611; -0.997248749443286], 1e-7);

%!test
%! ## a p plus one column terms matrix produces a model with an intercept
%! m = fitlm (X, y, [0 0 0; 1 0 0; 0 1 0]);
%! assert (m.NumCoefficients, 3);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.Coefficients.Estimate, ...
%!         [0.116188677790207; 2.50845149057086; -0.978835329825186], 1e-7);

%!test
%! ## a table with a Wilkinson formula fits the same model and predicts on a table
%! T = table (X(:,1), X(:,2), y, 'VariableNames', {'a','b','resp'});
%! m = fitlm (T, 'resp ~ a + b');
%! assert (m.NumCoefficients, 3);
%! assert (m.ResponseName, 'resp');
%! assert (m.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-8);
%! Xt = table ([0.5;1.0], [0.25;1.0], 'VariableNames', {'a','b'});
%! yp = predict (m, Xt);
%! assert (yp(1), 1.125705590619342, 1e-10);
%! assert (yp(2), 1.645804838535884, 1e-10);

%!test
%! ## a matrix with a Wilkinson formula string fits the same model as the matrix alone
%! m = fitlm (X, y, 'y ~ x1 + x2');
%! assert (m.NumCoefficients, 3);
%! assert (m.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-8);

%!test
%! ## a table input with the default formula fits the same model as the matrix
%! T3 = table (X(:,1), X(:,2), y, 'VariableNames', {'x1','x2','y'});
%! m = fitlm (T3);
%! assert (m.ResponseName, 'y');
%! assert (m.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-8);

%!test
%! ## VarNames sets custom names and ResponseVar overrides the response name
%! m  = fitlm (X, y, 'VarNames', {'alpha','beta','resp'});
%! assert (m.ResponseName, 'resp');
%! assert (isequal (m.PredictorNames, {'alpha','beta'}));
%! assert (any (strcmp (m.CoefficientNames, 'alpha')));
%! assert (any (strcmp (m.CoefficientNames, 'beta')));
%! m2 = fitlm (X, y, 'VarNames', {'a','b','r'}, 'ResponseVar', 'r');
%! assert (m2.ResponseName, 'r');

%!test
%! ## a rank deficient design matrix leaves the dropped coefficients as NaN across the board
%! X_rd = [ones(n,1), X, X(:,1)+X(:,2)];
%! m = fitlm (X_rd, y);
%! assert (m.NumCoefficients, 5);
%! assert (m.NumEstimatedCoefficients, 3);
%! drop = find (m.Coefficients.SE == 0);
%! assert (numel (drop), 2);
%! assert (all (isnan (m.Coefficients.tStat(drop))));
%! assert (all (isnan (m.Coefficients.pValue(drop))));
%! assert (m.SST, 5.839104200023459e+02, 1e-8);
%! assert (all (all (m.CoefficientCovariance(drop,:) == 0)));
%! yp = predict (m, X_rd);
%! assert (size (yp), [n, 1]);
%! assert (! any (isnan (yp)));
%! assert (size (m.Diagnostics.Dfbetas), [20, 5]);
%! assert (all (isnan (m.Diagnostics.Dfbetas(:, drop)(:))));
%! assert (m.Diagnostics.Leverage(1), 0.370779220779221, 1e-8);

%!test
%! ## Intercept=false
%! mni = fitlm (X, y, 'Intercept', false);
%! assert (mni.NumCoefficients, 2);
%! assert (mni.Formula.HasIntercept, false);
%! assert (! any (strcmp (mni.CoefficientNames, '(Intercept)')));
%! [yp, yci] = predict (mni, [0.5 0.25; 1.0 1.0]);
%! assert (yp(1), 1.231398619227234, 1e-10);
%! assert (yp(2), 1.964172863732825, 1e-10);
%! assert (yci(1,1), 1.001262470857215, 1e-10);

%!test
%! ## p-column terms matrix
%! m_p = fitlm (X, y, [1 0; 0 1]);
%! assert (m_p.NumCoefficients, 2);
%! assert (! any (strcmp (m_p.CoefficientNames, '(Intercept)')));

%!test
%! ## p+1 column terms matrix
%! m_p1 = fitlm (X, y, [0 0 0; 1 0 0; 0 1 0]);
%! assert (m_p1.NumCoefficients, 3);
%! assert (m_p1.CoefficientNames{1}, '(Intercept)');

%!test
%! ## table Wilkinson formula
%! T = table (X(:,1), X(:,2), y, 'VariableNames', {'a','b','resp'});
%! mf = fitlm (T, 'resp ~ a + b');
%! assert (mf.NumCoefficients, 3);
%! assert (mf.ResponseName, 'resp');
%! assert (mf.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-8);
%! Xt = table ([0.5;1.0], [0.25;1.0], 'VariableNames', {'a','b'});
%! yp = predict (mf, Xt);
%! assert (yp(1), 1.125705590619342, 1e-10);
%! assert (yp(2), 1.645804838535884, 1e-10);

%!test
%! ## matrix Wilkinson formula
%! mfm = fitlm (X, y, 'y ~ x1 + x2');
%! assert (mfm.NumCoefficients, 3);
%! assert (mfm.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-8);

%!test
%! ## table default
%! T3 = table (X(:,1), X(:,2), y, 'VariableNames', {'x1','x2','y'});
%! mt = fitlm (T3);
%! assert (mt.ResponseName, 'y');
%! assert (mt.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-8);

%!test
%! ## VarNames sets custom names
%! vn = fitlm (X, y, 'VarNames', {'alpha','beta','resp'});
%! assert (vn.ResponseName, 'resp');
%! assert (isequal (vn.PredictorNames, {'alpha','beta'}));
%! assert (any (strcmp (vn.CoefficientNames, 'alpha')));
%! assert (any (strcmp (vn.CoefficientNames, 'beta')));

%!test
%! ## ResponseVar overrides VarNames
%! rv = fitlm (X, y, 'VarNames', {'a','b','r'}, 'ResponseVar', 'r');
%! assert (rv.ResponseName, 'r');

%!test
%! ## rank-deficient matrix
%! X_rd = [ones(n,1), X, X(:,1)+X(:,2)];
%! m_rd = fitlm (X_rd, y);
%! assert (m_rd.NumCoefficients, 5);
%! assert (m_rd.NumEstimatedCoefficients, 3);
%! drop = find (m_rd.Coefficients.SE == 0);
%! assert (numel (drop), 2);
%! assert (all (isnan (m_rd.Coefficients.tStat(drop))));
%! assert (all (isnan (m_rd.Coefficients.pValue(drop))));
%! assert (m_rd.SST, 5.839104200023459e+02, 1e-8);
%! assert (all (all (m_rd.CoefficientCovariance(drop,:) == 0)));
%! yp = predict (m_rd, X_rd);
%! assert (size (yp), [n, 1]);
%! assert (! any (isnan (yp)));
%! assert (yp(1:5), [0.192669485827486; 0.171266760882252; ...
%!                   0.0519805029545;  -0.165189287955771; ...
%!                  -0.480242611848561], 1e-10);

%!test
%! ## predict: ypred and default CI at new points
%! [yp, yci] = predict (mdl, [0.5 0.25; 1.0 1.0]);
%! assert (yp(1),    1.125705590619347, 1e-10);
%! assert (yp(2),    1.645804838535894, 1e-10);
%! assert (yci(1,1), 0.810180780547215, 1e-10);
%! assert (yci(1,2), 1.441230400691478, 1e-10);
%! assert (yci(2,1), 0.858229321851723, 1e-10);
%! assert (yci(2,2), 2.433380355220066, 1e-10);

%!test
%! ## predict: observation interval
%! [~, yci] = predict (mdl, [0.5 0.25; 1.0 1.0], 'Prediction', 'observation');
%! assert (yci(1,1), 0.677632064105988, 1e-10);
%! assert (yci(1,2), 1.573779117132706, 1e-10);

%!test
%! ## predict: alpha 0.01
%! [~, yci] = predict (mdl, [0.5 0.25; 1.0 1.0], 'Alpha', 0.01);
%! assert (yci(1,1), 0.692272619570008, 1e-10);
%! assert (yci(1,2), 1.559138561668685, 1e-10);

%!test
%! ## predict: simultaneous CI
%! [~, yci] = predict (mdl, [0.5 0.25; 1.0 1.0], 'Simultaneous', true);
%! assert (yci(1,1), 0.662572505689338, 1e-10);
%! assert (yci(1,2), 1.588838675549355, 1e-10);

%!test
%! ## predict: no Xnew returns all rows including training
%! [yp, yci] = predict (mdl);
%! assert (size (yp),  [20, 1]);
%! assert (size (yci), [20, 2]);
%! assert (yp(1),    0.192669485827490, 1e-10);
%! assert (yp(2),    0.171266760882255, 1e-10);
%! assert (yci(1,1), -0.001052067982566, 1e-10);
%! assert (yci(1,2),  0.386391039637546, 1e-10);

%!test
%! ## predict: NaN predictor propagates to NaN output and CI
%! [yp, yci] = predict (mdl, [0.5 0.25; NaN 1.0; 1.0 1.0]);
%! assert (yp(1), 1.125705590619347, 1e-10);
%! assert (isnan (yp(2)));
%! assert (yp(3), 1.645804838535894, 1e-10);
%! assert (isnan (yci(2,1)));
%! assert (isnan (yci(2,2)));

%!test
%! ## predict: categorical model predictions at group centres
%! Xc    = [1;1;1;2;2;2;3;3;3];
%! yc    = [2.1;2.3;1.9; 4.1;3.9;4.2; 6.3;5.8;6.1];
%! m_cat = fitlm (Xc, yc, 'linear', 'CategoricalVars', 1);
%! [yp, yci] = predict (m_cat, [1;2;3]);
%! assert (yp(1), 2.099999999999998, 1e-10);
%! assert (yp(2), 4.066666666666667, 1e-10);
%! assert (yp(3), 6.066666666666666, 1e-10);
%! assert (yci(1,1), 1.80971256321669, 1e-10);
%! assert (yci(1,2), 2.3902874367833,  1e-10);
%! assert (yci(2,1), 3.77637922988336, 1e-10);
%! assert (yci(2,2), 4.35695410344997, 1e-10);
%! assert (yci(3,1), 5.77637922988336, 1e-10);
%! assert (yci(3,2), 6.35695410344997, 1e-10);

%!test
%! ## predict: interaction model
%! [yp, yci] = predict (fitlm (X, y, 'interactions'), [0.5 0.25; 1.0 1.0]);
%! assert (yp(1),    0.964032452046850, 1e-10);
%! assert (yp(2),    1.282176811827644, 1e-10);
%! assert (yci(1,1), -0.110763003580605, 1e-10);
%! assert (yci(1,2),  2.038827907674306, 1e-10);

%!test
%! ## predict: weighted model, ypred and CI
%! w  = (1:n)' / sum (1:n);
%! mw = fitlm (X, y, 'Weights', w);
%! [yp, yci] = predict (mw, [0.5 0.25; 1.0 1.0]);
%! assert (yp(1),    1.15833357370544, 1e-10);
%! assert (yp(2),    1.74408669002694, 1e-10);
%! assert (yci(1,1), 0.802165170771357, 1e-10);
%! assert (yci(1,2), 1.51450197663953,  1e-10);
%! assert (yci(2,1), 0.69968979253134,  1e-10);
%! assert (yci(2,2), 2.78848358752254,  1e-10);

%!test
%! ## predict: no-intercept model, ypred and CI
%! mni = fitlm (X, y, 'Intercept', false);
%! [yp, yci] = predict (mni, [0.5 0.25; 1.0 1.0]);
%! assert (yp(1),    1.23139861922723, 1e-10);
%! assert (yp(2),    1.96417286373283, 1e-10);
%! assert (yci(1,1), 1.00126247085704, 1e-10);
%! assert (yci(1,2), 1.46153476759743, 1e-10);
%! assert (yci(2,1), 1.51833851162232, 1e-10);
%! assert (yci(2,2), 2.41000721584333, 1e-10);

%!test
%! ## predict: observation interval combined with simultaneous bound
%! [~, yci] = predict (mdl, [0.5 0.25; 1.0 1.0], ...
%!                      'Prediction', 'observation', 'Simultaneous', true);
%! assert (yci(1,1), 0.46801507632267,  1e-10);
%! assert (yci(1,2), 1.78339610491601,  1e-10);
%! assert (yci(2,1), 0.399032373599106, 1e-10);
%! assert (yci(2,2), 2.89257730347266,  1e-10);

%!test
%! ## output is 2x1 double column vector
%! ysim = random (mdl, [0.5, 0.25; 1.0, 1.0]);
%! assert (size (ysim), [2, 1]);
%! assert (class (ysim), 'double');
%! assert (iscolumn (ysim));

%!test
%! ## single row input gives 1x1 output
%! assert (size (random (mdl, [0.5, 0.25])), [1, 1]);

%!test
%! ## predict values are exact and noise added is finite
%! ypred = predict (mdl, [0.5, 0.25; 1.0, 1.0]);
%! ysim  = random (mdl, [0.5, 0.25; 1.0, 1.0]);
%! assert (ypred(1), 1.125705590619342, 1e-10);
%! assert (ypred(2), 1.645804838535884, 1e-10);
%! assert (all (isfinite (ysim - ypred)));

%!test
%! ## NaN predictor row gives NaN output, other rows stay finite
%! ysim = random (mdl, [0.5, 0.25; NaN, 1.0; 1.0, 1.0]);
%! assert (size (ysim), [3, 1]);
%! assert (isfinite (ysim(1)));
%! assert (isnan (ysim(2)));
%! assert (isfinite (ysim(3)));

%!test
%! ## two sequential calls produce different output
%! ya = random (mdl, [0.5, 0.25]);
%! yb = random (mdl, [0.5, 0.25]);
%! assert (! isequal (ya, yb));

%!test
%! ## random: table input, full training data, weighted and no-intercept
%! ## models all give finite output of the expected size
%! Xt = table ([0.5;1.0], [0.25;1.0], 'VariableNames', {'x1','x2'});
%! mw  = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! mni = fitlm (X, y, 'Intercept', false);
%! assert (size (random (mdl, Xt)), [2, 1]);
%! assert (all (isfinite (random (mdl, Xt))));
%! assert (size (random (mdl, X)), [20, 1]);
%! assert (sum (isnan (random (mdl, X))), 0);
%! assert (all (isfinite (random (mw,  [0.5, 0.25; 1.0, 1.0]))));
%! assert (all (isfinite (random (mni, [0.5, 0.25; 1.0, 1.0]))));

%!test
%! yf = feval (mdl, [0.5 0.25; 1.0 1.0; 0.2 0.04]);
%! assert (size (yf), [3, 1]);
%! assert (class (yf), 'double');
%! assert (yf(1), 1.125705590619342, 1e-10);
%! assert (yf(2), 1.645804838535884, 1e-10);
%! assert (yf(3), 0.578725562711373, 1e-10);
%! assert (yf, predict (mdl, [0.5 0.25; 1.0 1.0; 0.2 0.04]), 1e-10);

%!test
%! yf = feval (mdl, [0.5; 1.0; 0.2], [0.25; 1.0; 0.04]);
%! assert (size (yf), [3, 1]);
%! assert (iscolumn (yf));
%! assert (yf, predict (mdl, [0.5 0.25; 1.0 1.0; 0.2 0.04]), 1e-10);

%!test
%! yf = feval (mdl, [0.5, 1.0, 0.2], [0.25, 1.0, 0.04]);
%! assert (size (yf), [1, 3]);
%! assert (isrow (yf));
%! assert (yf(:), predict (mdl, [0.5 0.25; 1.0 1.0; 0.2 0.04]), 1e-10);

%!test
%! yf = feval (mdl, 0.5, 0.25);
%! assert (size (yf), [1, 1]);
%! assert (yf, 1.125705590619342, 1e-10);
%! assert (yf, predict (mdl, [0.5 0.25]), 1e-10);

%!test
%! yf = feval (mdl, 0.5, [0.1; 0.2; 0.3]);
%! assert (size (yf), [3, 1]);
%! assert (yf(1), 1.272530890093120, 1e-10);
%! assert (yf(2), 1.174647357110602, 1e-10);
%! assert (yf(3), 1.076763824128083, 1e-10);
%! assert (yf, predict (mdl, [0.5 0.1; 0.5 0.2; 0.5 0.3]), 1e-10);

%!test
%! yf = feval (mdl, [0.1; 0.5; 0.9], 0.25);
%! assert (size (yf), [3, 1]);
%! assert (yf(1), 0.122324994390997, 1e-10);
%! assert (yf(2), 1.125705590619342, 1e-10);
%! assert (yf(3), 2.129086186847688, 1e-10);
%! assert (yf, predict (mdl, [0.1 0.25; 0.5 0.25; 0.9 0.25]), 1e-10);

%!test
%! m = fitlm ((1:n)' / n, 2 * (1:n)' / n + 0.1 * sin ((1:n)'));
%! assert (size (feval (m, 0.5)), [1, 1]);
%! assert (size (feval (m, [0.3; 0.5; 0.9])), [3, 1]);
%! assert (feval (m, 0.5), predict (m, 0.5), 1e-10);
%! assert (feval (m, [0.3; 0.5; 0.9]), predict (m, [0.3; 0.5; 0.9]), 1e-10);

%!test
%! T  = table ([0.5; 1.0; 0.2], [0.25; 1.0; 0.04], 'VariableNames', {'x1', 'x2'});
%! yf = feval (mdl, T);
%! assert (size (yf), [3, 1]);
%! assert (yf, predict (mdl, [0.5 0.25; 1.0 1.0; 0.2 0.04]), 1e-10);

%!test
%! yf = feval (mdl, [0.5 0.25; NaN 1.0; 1.0 1.0]);
%! assert (isfinite (yf(1)));
%! assert (isnan (yf(2)));
%! assert (isfinite (yf(3)));

%!test
%! yf = feval (mdl, [0.5; NaN; 1.0], [0.25; 1.0; 1.0]);
%! assert (isnan (yf(2)));
%! yf = feval (mdl, [0.5; 1.0; 1.0], [0.25; NaN; 1.0]);
%! assert (isnan (yf(2)));

%!test
%! yf = feval (mdl, X);
%! assert (size (yf), [20, 1]);
%! assert (yf, mdl.Fitted, 1e-10);

%!test
%! m  = fitlm (X, y, 'Intercept', false);
%! yf = feval (m, [0.5 0.25; 1.0 1.0]);
%! assert (yf, predict (m, [0.5 0.25; 1.0 1.0]), 1e-10);
%! assert (feval (m, [0.5; 1.0], [0.25; 1.0]), yf, 1e-10);

%!test
%! m  = fitlm (X, y, 'interactions');
%! yf = feval (m, [0.5 0.25; 1.0 1.0]);
%! assert (yf, predict (m, [0.5 0.25; 1.0 1.0]), 1e-10);
%! assert (feval (m, [0.5; 1.0], [0.25; 1.0]), yf, 1e-10);

%!test
%! m  = fitlm ([1;1;1;2;2;2;3;3;3], [2.1;2.3;1.9;4.1;3.9;4.2;6.3;5.8;6.1], ...
%!             'linear', 'CategoricalVars', 1);
%! yf = feval (m, [1; 2; 3]);
%! assert (yf(1), 2.099999999999998, 1e-10);
%! assert (yf(2), 4.066666666666667, 1e-10);
%! assert (yf(3), 6.066666666666666, 1e-10);

%!test
%! ci = coefCI (mdl);
%! assert (size (ci), [3, 2]);
%! assert (class (ci), 'double');
%! assert (all (ci(:,1) < ci(:,2)));
%! assert (ci(1,1), -0.120502736154050,  1e-10);
%! assert (ci(1,2),  0.352880091734465, 1e-10);
%! assert (ci(2,1),  1.470249604061007,  1e-10);
%! assert (ci(2,2),  3.546653377080718,  1e-10);
%! assert (ci(3,1), -1.026857022014626,  1e-10);
%! assert (ci(3,2), -0.930813637635746, 1e-10);

%!test
%! ## midpoints equal estimates
%! ci = coefCI (mdl);
%! t  = tinv (0.975, mdl.DFE);
%! assert ((ci(:,1) + ci(:,2)) / 2, mdl.Coefficients.Estimate, 1e-10);
%! assert (ci(:,2) - ci(:,1), 2 * t * mdl.Coefficients.SE, 1e-10);

%!test
%! assert (coefCI (mdl, 0.05), coefCI (mdl));

%!test
%! ci   = coefCI (mdl);
%! ci01 = coefCI (mdl, 0.01);
%! t01  = tinv (0.995, mdl.DFE);
%! assert (size (ci01), [3, 2]);
%! assert (ci01(1,1), -0.208951721610638, 1e-10);
%! assert (ci01(1,2),  0.441329077191052, 1e-10);
%! assert (ci01(2,1),  1.08228494564489,  1e-10);
%! assert (ci01(2,2),  3.934618035496833, 1e-10);
%! assert (ci01(3,1), -1.044802201703589, 1e-10);
%! assert (ci01(3,2), -0.912868457946783, 1e-10);
%! assert (all ((ci01(:,2) - ci01(:,1)) > (ci(:,2) - ci(:,1))));
%! assert (ci01(:,2) - ci01(:,1), 2 * t01 * mdl.Coefficients.SE, 1e-10);

%!test
%! ci0 = coefCI (mdl, 0);
%! assert (all (ci0(:,1) == -Inf));
%! assert (all (ci0(:,2) == +Inf));

%!test
%! ## alpha=1 collapses to point estimates
%! ci1 = coefCI (mdl, 1);
%! assert (ci1(:,1), mdl.Coefficients.Estimate, 1e-10);
%! assert (ci1(:,2), mdl.Coefficients.Estimate, 1e-10);

%!test
%! m  = fitlm (X, y, 'Intercept', false);
%! ci = coefCI (m);
%! t  = tinv (0.975, m.DFE);
%! assert (size (ci), [2, 2]);
%! assert (ci(1,1), 2.486679110991696,  1e-10);
%! assert (ci(1,2), 3.436164115360526, 1e-10);
%! assert (ci(2,1), -1.027166590567854, 1e-10);
%! assert (ci(2,2), -0.967330908318718, 1e-10);
%! assert ((ci(:,1) + ci(:,2)) / 2, m.Coefficients.Estimate, 1e-10);
%! assert (ci(:,2) - ci(:,1), 2 * t * m.Coefficients.SE, 1e-10);

%!test
%! m  = fitlm (X, y, 'interactions');
%! ci = coefCI (m);
%! t  = tinv (0.975, m.DFE);
%! assert (size (ci), [4, 2]);
%! assert (ci(1,1), -0.201030907566802, 1e-10);
%! assert (ci(1,2),  0.516312363642881, 1e-10);
%! assert ((ci(:,1) + ci(:,2)) / 2, m.Coefficients.Estimate, 1e-10);
%! assert (ci(:,2) - ci(:,1), 2 * t * m.Coefficients.SE, 1e-10);

%!test
%! ## constant model (1 coefficient)
%! m  = fitlm (X, y, 'constant');
%! ci = coefCI (m);
%! t  = tinv (0.975, m.DFE);
%! assert (size (ci), [1, 2]);
%! assert ((ci(1,1) + ci(1,2)) / 2, m.Coefficients.Estimate, 1e-10);
%! assert (ci(1,2) - ci(1,1), 2 * t * m.Coefficients.SE, 1e-10);

%!test
%! m  = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! ci = coefCI (m);
%! t  = tinv (0.975, m.DFE);
%! assert (size (ci), [3, 2]);
%! assert (ci(1,1), -0.355978167660141, 1e-10);
%! assert (ci(1,2),  0.516619434992026, 1e-10);
%! assert ((ci(:,1) + ci(:,2)) / 2, m.Coefficients.Estimate, 1e-10);
%! assert (ci(:,2) - ci(:,1), 2 * t * m.Coefficients.SE, 1e-10);

%!test
%! ## rank-deficient: dropped rows give [0,0], active rows are finite
%! m    = fitlm ([ones(n,1), X, X(:,1)+X(:,2)], y);
%! ci   = coefCI (m);
%! drop = find (m.Coefficients.SE == 0);
%! assert (size (ci), [5, 2]);
%! assert (all (all (ci(drop, :) == 0)));
%! assert (all (all (isfinite (ci(setdiff (1:5, drop'), :)))));

%!test
%! m  = fitlm ([1;1;1;2;2;2;3;3;3], [2.1;2.3;1.9;4.1;3.9;4.2;6.3;5.8;6.1], ...
%!             'linear', 'CategoricalVars', 1);
%! ci = coefCI (m);
%! assert (size (ci), [3, 2]);
%! assert (ci(1,1), 1.80971256321669,  1e-10);
%! assert (ci(1,2), 2.3902874367833,   1e-10);
%! assert (ci(2,1), 1.55613823658119,  1e-10);
%! assert (ci(2,2), 2.37719509675214,  1e-10);
%! assert (ci(3,1), 3.55613823658119,  1e-10);
%! assert (ci(3,2), 4.37719509675214,  1e-10);

%!test
%! [p, F, r] = coefTest (mdl);
%! assert (size (p), [1, 1]);
%! assert (class (p), 'double');
%! assert (p >= 0 && p <= 1);
%! assert (F >= 0);
%! assert (p, 9.489880832170599e-28, -1e-8);
%! assert (F, 1.283149098426142e+04, -1e-8);
%! assert (r, 2);

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
%! assert (F, Fm, -1e-10);
%! assert (p, pm, -1e-10);

%!test
%! ## explicit H matches default
%! k     = mdl.NumCoefficients;
%! H_exp = [zeros(k-1, 1), eye(k-1)];
%! [p1, F1, r1] = coefTest (mdl);
%! [p2, F2, r2] = coefTest (mdl, H_exp);
%! assert (p2, p1, -1e-10);
%! assert (F2, F1, -1e-10);
%! assert (r2, size (H_exp, 1));

%!test
%! ## pinned single and joint H
%! [p1, F1, r1] = coefTest (mdl, [1 0 0]);
%! assert (p1, 0.314859866747774, -1e-8);
%! assert (F1, 1.072634101844537, -1e-8);
%! assert (r1, 1);
%! [p2, F2, r2] = coefTest (mdl, [0 1 0]);
%! assert (p2, 8.937794169018252e-05, -1e-8);
%! assert (F2, 25.985840929474932, -1e-8);
%! assert (r2, 1);
%! [p3, F3, r3] = coefTest (mdl, [0 0 1]);
%! assert (p3, 8.656938305821102e-19, -1e-8);
%! assert (F3, 1.849410599855684e+03, -1e-8);
%! assert (r3, 1);
%! [pm, Fm, rm] = coefTest (mdl, [0 1 0; 0 0 1]);
%! assert (pm, 9.489880832170599e-28, -1e-8);
%! assert (Fm, 1.283149098426142e+04, -1e-8);
%! assert (rm, 2);

%!test
%! ## trivial hypothesis and C=0
%! b        = mdl.Coefficients.Estimate;
%! [p0, F0] = coefTest (mdl, [0 1 0], b(2));
%! assert (F0 < 1e-12);
%! assert (p0, 1, 1e-10);
%! [pa, Fa] = coefTest (mdl, [0 1 0], 0);
%! [pb, Fb] = coefTest (mdl, [0 1 0]);
%! assert (pa, pb, -1e-10);
%! assert (Fa, Fb, -1e-10);

%!test
%! ## H with C
%! [pc, Fc, rc] = coefTest (mdl, [0 1 0; 0 0 1], [1.5; -1.0]);
%! assert (pc, 2.833788304242915e-09, -1e-8);
%! assert (Fc, 77.603887650386312, -1e-8);
%! assert (rc, 2);
%! [pr, Fr] = coefTest (mdl, [0 1 0; 0 0 1], [1.5, -1.0]);
%! assert (pr, pc, -1e-10);
%! assert (Fr, Fc, -1e-10);
%! [ps, Fs] = coefTest (mdl, [0 1 0], 1.5);
%! assert (ps, 0.056184159363707, -1e-8);
%! assert (Fs, 4.199865537706047, -1e-8);

%!test
%! ## no-intercept model
%! m = fitlm (X, y, 'Intercept', false);
%! [p, F, r] = coefTest (m);
%! assert (p, 6.060655830723051e-32, -1e-8);
%! assert (F, 2.646694317541346e+04, -1e-8);
%! assert (r, m.NumCoefficients);
%! [p2, F2] = coefTest (m, eye (m.NumCoefficients));
%! assert (p2, p, -1e-10);
%! assert (F2, F, -1e-10);

%!test
%! ## interaction model
%! m = fitlm (X, y, 'interactions');
%! [p, F, r] = coefTest (m);
%! assert (p, 1.164196605688161e-25, -1e-8);
%! assert (F, 8.107508574885546e+03, -1e-8);
%! assert (r, m.NumCoefficients - 1);
%! assert (r != m.NumPredictors);

%!test
%! ## weighted model
%! m = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! [p, F, r] = coefTest (m);
%! assert (p, 1.481920976389473e-27, -1e-8);
%! assert (F, 1.217557180481257e+04, -1e-8);
%! assert (r, 2);
%! assert (p, m.ModelFitVsNullModel.Pvalue, -1e-8);

%!test
%! ## categorical model
%! m = fitlm ([1;1;1;2;2;2;3;3;3], [2.1;2.3;1.9;4.1;3.9;4.2;6.3;5.8;6.1], ...
%!            'linear', 'CategoricalVars', 1);
%! [p, F, r] = coefTest (m);
%! assert (p, 1.197590680415813e-06, -1e-8);
%! assert (F, 2.795000000000035e+02, -1e-8);
%! assert (r, 2);
%! [p1, F1] = coefTest (m, [1 0 0]);
%! assert (F1, 3.133421052631613e+02, -1e-8);
%! assert (p1, 2.087464608380450e-06, -1e-8);
%! [p2, F2] = coefTest (m, [0 1 0]);
%! assert (F2, 1.374078947368438e+02, -1e-8);
%! assert (p2, 2.325514143662469e-05, -1e-8);
%! [p3, F3] = coefTest (m, [0 0 1]);
%! assert (F3, 5.589868421052698e+02, -1e-8);
%! assert (p3, 3.757733067786492e-07, -1e-8);

%!test
%! ## constant model
%! m = fitlm (X, y, 'constant');
%! [p, F, r] = coefTest (m);
%! assert (p, 0.000239936408695073, -1e-8);
%! assert (F, 20.3359164947506, -1e-8);
%! assert (r, 1);

%!test
%! ## rank-deficient model
%! m    = fitlm ([ones(n,1), X, X(:,1)+X(:,2)], y);
%! [p, F] = coefTest (m);
%! assert (isnan (p));
%! assert (isnan (F));
%! drop = find (m.Coefficients.SE == 0);
%! keep = setdiff (2:m.NumCoefficients, drop');
%! H    = zeros (numel (keep), m.NumCoefficients);
%! for i = 1:numel (keep)
%!   H(i, keep(i)) = 1;
%! endfor
%! [p2, F2, r2] = coefTest (m, H);
%! assert (p2, 6.70657058643085e-30, -1e-8);
%! assert (F2, 17716.1864263456, -1e-8);
%! assert (r2, numel (keep));

%!test
%! p = dwtest (mdl);
%! assert (size (p),  [1, 1]);
%! assert (class (p), 'double');
%! [p, DW] = dwtest (mdl);
%! assert (size (DW), [1, 1]);
%! assert (p  >= 0 && p  <= 1);
%! assert (DW >= 0 && DW <= 4);
%! assert (p,  4.702593821571290e-04, -1e-6);
%! assert (DW, 0.870000704251173, 1e-12);

%!test
%! [p1, DW1] = dwtest (mdl);
%! [p2, DW2] = dwtest (mdl, 'exact', 'both');
%! assert (p1, p2, 1e-14);
%! assert (DW1, DW2, 1e-14);

%!test
%! ## DW is the same for all method and tail options
%! [~, d1] = dwtest (mdl, 'exact', 'both');
%! [~, d2] = dwtest (mdl, 'exact', 'right');
%! [~, d3] = dwtest (mdl, 'exact', 'left');
%! [~, d4] = dwtest (mdl, 'approximate', 'both');
%! [~, d5] = dwtest (mdl, 'approximate', 'right');
%! [~, d6] = dwtest (mdl, 'approximate', 'left');
%! assert (d1, 0.870000704251173, 1e-12);
%! assert (d2, 0.870000704251173, 1e-12);
%! assert (d3, 0.870000704251173, 1e-12);
%! assert (d4, 0.870000704251173, 1e-12);
%! assert (d5, 0.870000704251173, 1e-12);
%! assert (d6, 0.870000704251173, 1e-12);

%!test
%! ## one-sided p-values sum to 1 and two-sided equals twice the smaller
%! pb = dwtest (mdl, 'exact', 'both');
%! pr = dwtest (mdl, 'exact', 'right');
%! pl = dwtest (mdl, 'exact', 'left');
%! assert (pr + pl, 1, 1e-12);
%! assert (pb, 4.702593821571290e-04, 1e-12);

%!test
%! ## all six method and tail combinations pinned
%! assert (dwtest (mdl, 'exact', 'both'), 4.702593821571290e-04, -1e-6);
%! assert (dwtest (mdl, 'exact', 'right'), 2.351296910785645e-04, -1e-6);
%! assert (dwtest (mdl, 'exact', 'left'), 0.999764870308921, -1e-6);
%! assert (dwtest (mdl, 'approximate', 'both'), 0.001058795514879, -1e-6);
%! assert (dwtest (mdl, 'approximate', 'right'), 5.293977574395035e-04, -1e-6);
%! assert (dwtest (mdl, 'approximate', 'left'), 0.999470602242560, -1e-6);

%!test
%! ## no-intercept model
%! m = fitlm (X, y, 'Intercept', false);
%! [p, DW] = dwtest (m, 'exact', 'both');
%! assert (DW, 0.841468411374128, 1e-12);
%! assert (p, 0.001402191159200, -1e-6);
%! assert (dwtest (m, 'exact', 'right'), 7.010955795999754e-04, -1e-6);
%! assert (dwtest (m, 'approximate', 'right'), 0.001350534002321, -1e-6);

%!test
%! ## weighted model
%! m = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! [p, DW] = dwtest (m, 'exact', 'both');
%! assert (DW, 0.871162354803032, 1e-12);
%! assert (p, 4.771641146603785e-04, -1e-6);
%! assert (dwtest (m, 'exact', 'right'), 2.385820573301892e-04, -1e-6);
%! assert (dwtest (m, 'approximate', 'right'), 5.346779629058873e-04, -1e-6);

%!test
%! ## positive autocorrelation model
%! m = fitlm ((1:n)'/n, sin (pi * (1:n)'/n));
%! [~, DW] = dwtest (m, 'exact', 'both');
%! pr = dwtest (m, 'exact', 'right');
%! pl = dwtest (m, 'exact', 'left');
%! assert (DW, 0.118112272685229, 1e-10);
%! assert (DW < 1);
%! assert (pr < pl);
%! assert (pr < 1e-10);

%!test
%! ## negative autocorrelation model
%! m = fitlm ((1:n)'/n, 2*(1:n)'/n + repmat ([1; -1], n/2, 1));
%! [pb, DW] = dwtest (m, 'exact', 'both');
%! pl = dwtest (m, 'exact', 'left');
%! pr = dwtest (m, 'exact', 'right');
%! assert (pb, 4.205713999283489e-09, 1e-10);
%! assert (DW, 3.825974025974026, 1e-10);
%! assert (DW > 2);
%! assert (pl < pr);
%! assert (pb, 2 * pl, 1e-10);
%! assert (pb < 1e-7);

%!test
%! m = addTerms (mdl, 'x1:x2');
%! assert (isa (m, 'LinearModel'));
%! assert (mdl.NumCoefficients, 3);
%! assert (m.NumCoefficients, 4);
%! assert (m.NumPredictors, 2);
%! assert (m.NumObservations, 20);
%! assert (m.DFE, 16);
%! assert (m.Coefficients.Estimate(1), 0.157640728038039, -1e-8);
%! assert (m.Coefficients.Estimate(2), 2.085426803117909, -1e-8);
%! assert (m.Coefficients.Estimate(3), -0.929682701072813, -1e-8);
%! assert (m.Coefficients.Estimate(4), -0.031208018255475, -1e-8);
%! assert (m.Coefficients.SE(1), 0.169192291625763, -1e-8);
%! assert (m.Coefficients.SE(2), 1.361534257888685, -1e-8);
%! assert (m.Coefficients.SE(3), 0.148744319911833, -1e-8);
%! assert (m.Coefficients.SE(4), 0.0932669056882381, -1e-8);
%! assert (m.Coefficients.tStat(1), 0.931725237144526, -1e-8);
%! assert (m.Coefficients.tStat(2), 1.531674132351069, -1e-8);
%! assert (m.Coefficients.tStat(3), -6.250206405353000, -1e-8);
%! assert (m.Coefficients.tStat(4), -0.334609774230031, -1e-8);
%! assert (m.Coefficients.pValue(1), 0.365325503492671, -1e-8);
%! assert (m.Coefficients.pValue(2), 0.145134783727025, -1e-8);
%! assert (m.Coefficients.pValue(3), 1.159217784590233e-05, -1e-8);
%! assert (m.Coefficients.pValue(4), 0.742265736761240, -1e-8);
%! assert (m.SSE, 0.383859187927621, -1e-8);
%! assert (m.MSE, 0.023991199245515, -1e-8);
%! assert (m.RMSE, 0.154890926930905, -1e-8);
%! assert (m.Rsquared.Ordinary, 0.999342606032059, -1e-8);
%! assert (m.Rsquared.Adjusted, 0.999219344663070, -1e-8);
%! assert (m.LogLikelihood, 11.153346988927943, -1e-8);
%! assert (m.ModelFitVsNullModel.Fstat, 8.107508574898859e+03, -1e-6);
%! assert (m.ModelFitVsNullModel.Pvalue, 1.164196605672873e-25, -1e-6);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.CoefficientNames{3}, 'x2');
%! assert (m.CoefficientNames{4}, 'x1:x2');

%!test
%! ## x1*x2 crossing gives same result as x1:x2 when main effects exist
%! m = addTerms (mdl, 'x1*x2');
%! assert (m.NumCoefficients, 4);
%! assert (m.DFE, 16);
%! assert (m.SSE, 0.383859187927621, -1e-8);
%! assert (m.Coefficients.Estimate(1), 0.157640728038039, -1e-8);
%! assert (m.Coefficients.Estimate(2), 2.085426803117909, -1e-8);
%! assert (m.Coefficients.Estimate(3), -0.929682701072813, -1e-8);
%! assert (m.Coefficients.Estimate(4), -0.031208018255475, -1e-8);
%! assert (m.CoefficientNames{4}, 'x1:x2');

%!test
%! m = addTerms (mdl, 'x1 + x1:x2');
%! assert (m.NumCoefficients, 4);
%! assert (m.DFE, 16);
%! assert (m.SSE, 0.383859187927621, -1e-8);
%! assert (m.Coefficients.Estimate(1), 0.157640728038039, -1e-8);
%! assert (m.Coefficients.Estimate(2), 2.085426803117909, -1e-8);
%! assert (m.Coefficients.Estimate(3), -0.929682701072813, -1e-8);
%! assert (m.Coefficients.Estimate(4), -0.031208018255475, -1e-8);

%!test
%! ## adding existing term returns equivalent model
%! ws = warning ('off', 'all');
%! m  = addTerms (mdl, 'x1');
%! warning (ws);
%! assert (m.NumCoefficients, 3);
%! assert (m.DFE, 17);
%! assert (m.Coefficients.Estimate(1), 0.116188677790207, 1e-7);
%! assert (m.Coefficients.Estimate(2), 2.508451490570863, 1e-7);
%! assert (m.Coefficients.Estimate(3), -0.978835329825186, 1e-7);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.CoefficientNames{3}, 'x2');

%!test
%! m = addTerms (mdl, 'x2^2');
%! assert (m.NumCoefficients, 4);
%! assert (m.DFE, 16);
%! assert (m.SSE, 0.386103933724971, -1e-8);
%! assert (m.Coefficients.Estimate(1), 0.130152473216993, -1e-8);
%! assert (m.Coefficients.Estimate(2), 2.380771990884563, -1e-8);
%! assert (m.Coefficients.Estimate(3), -0.967672823484773, -1e-8);
%! assert (m.Coefficients.Estimate(4), -2.991483322469043e-04, -1e-8);
%! assert (m.Coefficients.SE(1), 0.154974488176692, -1e-8);
%! assert (m.Coefficients.SE(2), 1.071554310049276, -1e-8);
%! assert (m.Coefficients.SE(3), 0.085801310569364, -1e-8);
%! assert (m.Coefficients.SE(4), 0.002211890858232, -1e-8);
%! assert (m.CoefficientNames{4}, 'x2^2');

%!test
%! ## x1^2 rank-deficient: DFE unchanged SE zero for dropped term
%! m = addTerms (mdl, 'x1^2');
%! assert (m.NumCoefficients, 4);
%! assert (m.DFE, 17);
%! assert (m.SSE, 0.386545331386823, -1e-8);
%! assert (m.Coefficients.Estimate(4), 0);
%! assert (m.Coefficients.SE(4), 0);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.CoefficientNames{3}, 'x2');
%! assert (m.CoefficientNames{4}, 'x1^2');


%!test
%! ## numeric matrix [1,1,0] same as string x1:x2
%! m = addTerms (mdl, [1, 1, 0]);
%! assert (m.NumCoefficients, 4);
%! assert (m.DFE, 16);
%! assert (m.SSE, 0.383859187927621, -1e-8);
%! assert (m.Coefficients.Estimate(1), 0.157640728038039, -1e-8);
%! assert (m.Coefficients.Estimate(2), 2.085426803117909, -1e-8);
%! assert (m.Coefficients.Estimate(3), -0.929682701072813, -1e-8);
%! assert (m.Coefficients.Estimate(4), -0.031208018255475, -1e-8);

%!test
%! ## numeric matrix [1,1] auto-padded to [1,1,0]
%! m = addTerms (mdl, [1, 1]);
%! assert (m.NumCoefficients, 4);
%! assert (m.DFE, 16);
%! assert (m.SSE, 0.383859187927621, -1e-8);
%! assert (m.Coefficients.Estimate(1), 0.157640728038039, -1e-8);
%! assert (m.Coefficients.Estimate(2), 2.085426803117909, -1e-8);
%! assert (m.Coefficients.Estimate(3), -0.929682701072813, -1e-8);
%! assert (m.Coefficients.Estimate(4), -0.031208018255475, -1e-8);

%!test
%! m = addTerms (mdl, [1, 1, 0; 0, 2, 0]);
%! assert (m.NumCoefficients, 5);
%! assert (m.DFE, 15);
%! assert (m.SSE, 0.315784637443501, -1e-8);
%! assert (m.CoefficientNames{4}, 'x1:x2');
%! assert (m.CoefficientNames{5}, 'x2^2');

%!test
%! mc = fitlm (X, y, 'constant');
%! m  = addTerms (mc, 'x1');
%! assert (m.NumCoefficients, 2);
%! assert (m.DFE, 18);
%! assert (m.Coefficients.Estimate(1), 3.884704697617172, -1e-8);
%! assert (m.Coefficients.Estimate(2), -18.047090435758047, -1e-8);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');

%!test
%! ## step from constant to full linear model
%! mc  = fitlm (X, y, 'constant');
%! mc1 = addTerms (mc, 'x1');
%! mc2 = addTerms (mc1, 'x2');
%! assert (mc2.NumCoefficients, 3);
%! assert (mc2.DFE, 17);
%! assert (mc2.Coefficients.Estimate(1), 0.116188677790207, 1e-7);
%! assert (mc2.Coefficients.Estimate(2), 2.508451490570863, 1e-7);
%! assert (mc2.Coefficients.Estimate(3), -0.978835329825186, 1e-7);
%! assert (mc2.CoefficientNames{1}, '(Intercept)');
%! assert (mc2.CoefficientNames{2}, 'x1');
%! assert (mc2.CoefficientNames{3}, 'x2');

%!test
%! ## adding intercept to no-intercept model
%! mni = fitlm (X, y, 'Intercept', false);
%! m   = addTerms (mni, '1');
%! assert (m.NumCoefficients, 3);
%! assert (m.DFE, 17);
%! assert (m.Coefficients.Estimate(1), 0.116188677790207, 1e-7);
%! assert (m.Coefficients.Estimate(2), 2.508451490570863, 1e-7);
%! assert (m.Coefficients.Estimate(3), -0.978835329825186, 1e-7);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.CoefficientNames{3}, 'x2');

%!test
%! ## weighted model weights preserved
%! mw = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! m  = addTerms (mw, 'x1:x2');
%! assert (m.NumCoefficients, 4);
%! assert (m.DFE, 16);
%! assert (m.SSE, 0.019230053719402, -1e-8);
%! assert (m.Coefficients.Estimate(1), -0.122645849510537, -1e-8);
%! assert (m.Coefficients.Estimate(2), 4.125799311652051, -1e-8);
%! assert (m.Coefficients.Estimate(3), -1.128467507844852, -1e-8);
%! assert (m.Coefficients.Estimate(4), 0.081921546574140, -1e-8);
%! assert (m.Coefficients.SE(1), 0.354926338845420, -1e-8);
%! assert (m.Coefficients.SE(2), 2.205951720932299, -1e-8);
%! assert (m.Coefficients.SE(3), 0.205033663966776, -1e-8);
%! assert (m.Coefficients.SE(4), 0.115523382309399, -1e-8);
%! assert (m.CoefficientNames{4}, 'x1:x2');

%!test
%! ## excluded observations preserved
%! me = fitlm (X, y, 'Exclude', [1, 2]);
%! m  = addTerms (me, 'x1:x2');
%! assert (m.NumObservations, 18);
%! assert (m.DFE, 14);
%! assert (m.NumCoefficients, 4);
%! assert (m.Coefficients.Estimate(1), -0.345521184099998, -1e-8);
%! assert (m.Coefficients.Estimate(2), 5.139185607268283, -1e-8);
%! assert (m.Coefficients.Estimate(3), -1.198851436671170, -1e-8);
%! assert (m.Coefficients.Estimate(4), 0.112619530301200, -1e-8);
%! assert (m.CoefficientNames{4}, 'x1:x2');

%!test
%! ## remove two predictors by string from a 4-predictor model
%! Xh = [7 26 6 60; 1 29 15 52; 11 56 8 20; 11 31 8 47; 7 52 6 33; ...
%!        11 55 9 22; 3 71 17 6; 1 31 22 44; 2 54 18 22; 21 47 4 26; ...
%!        1 40 23 34; 11 66 9 12; 10 68 8 12];
%! yh = [78.5;74.3;104.3;87.6;95.9;109.2;102.7;72.5;93.1;115.9;83.8;113.3;109.4];
%! m  = removeTerms (fitlm (Xh, yh), 'x3 + x4');
%! assert (m.NumCoefficients, 3);
%! assert (m.NumEstimatedCoefficients, 3);
%! assert (m.DFE, 10);
%! assert (m.NumObservations, 13);
%! assert (m.NumVariables, 5);
%! assert (m.Coefficients.Estimate(1), 52.577348882089481, -1e-8);
%! assert (m.Coefficients.Estimate(2), 1.468305742215555, -1e-8);
%! assert (m.Coefficients.Estimate(3), 0.662250491274645, -1e-8);
%! assert (m.Coefficients.SE(1), 2.286174334503340, -1e-8);
%! assert (m.Coefficients.SE(2), 0.121300923606266, -1e-8);
%! assert (m.Coefficients.SE(3), 0.045854721468522, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.Coefficients.tStat(1), 22.997961305305111, -1e-8);
%! assert (m.Coefficients.tStat(2), 12.104654264476748, -1e-8);
%! assert (m.Coefficients.tStat(3), 14.442362096327519, -1e-8);
%! assert (m.Coefficients.pValue(1), 5.456570901490983e-10, -1e-7);
%! assert (m.Coefficients.pValue(2), 2.692212179685427e-07, -1e-8);
%! assert (m.Coefficients.pValue(3), 5.028960315638413e-08, -1e-8);
%! assert (m.SSE, 57.904483176113658, -1e-8);
%! assert (m.RMSE, 2.40633503852047, -1e-8);
%! assert (m.MSE, 5.790448317611299, 1e-12);
%! assert (m.SST, 2.715763076923078e+03, 1e-8);
%! assert (m.Rsquared.Ordinary, 0.978678374535632, -1e-8);
%! assert (m.Rsquared.Adjusted, 0.974414049442758, -1e-8);
%! assert (size (m.CoefficientCovariance), [3, 3]);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.CoefficientNames{3}, 'x2');
%! assert (m.Formula.HasIntercept, true);
%! assert (m.Formula.LinearPredictor, '1 + x1 + x2');
%! assert (height (m.Diagnostics), 13);
%! assert (sum (m.Diagnostics.Leverage), 3, 1e-10);
%! assert (m.Residuals.Raw, yh - m.Fitted, 1e-10);

%!test
%! m = removeTerms (mdl, 'x2');
%! assert (m.NumCoefficients, 2);
%! assert (m.NumEstimatedCoefficients, 2);
%! assert (m.DFE, 18);
%! assert (m.NumObservations, 20);
%! assert (m.Coefficients.Estimate(1), 3.88470469761717, -1e-8);
%! assert (m.Coefficients.Estimate(2), -18.047090435758, -1e-8);
%! assert (m.Coefficients.SE(2), 1.19086428900602, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.Coefficients.tStat(2), -15.1546155194741, -1e-8);
%! assert (m.SSE, 42.4383708132815, -1e-8);
%! assert (m.SST, 583.910420002346, -1e-8);
%! assert (m.Rsquared.Ordinary, 0.927320408474452, -1e-8);
%! assert (size (m.CoefficientCovariance), [2, 2]);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.Formula.HasIntercept, true);
%! assert (m.Formula.LinearPredictor, '1 + x1');
%! assert (height (m.Diagnostics), 20);
%! assert (m.Residuals.Raw, y - m.Fitted, 1e-10);

%!test
%! ## removing the intercept via string '1'
%! m = removeTerms (mdl, '1');
%! assert (m.NumCoefficients, 2);
%! assert (m.NumEstimatedCoefficients, 2);
%! assert (m.DFE, 18);
%! assert (m.NumObservations, 20);
%! assert (m.Formula.HasIntercept, false);
%! assert (m.Formula.LinearPredictor, 'x1 + x2');
%! assert (m.Coefficients.Estimate(1), 2.96142161317611, -1e-8);
%! assert (m.Coefficients.Estimate(2), -0.997248749443286, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.SSE, 0.410934843407688, -1e-8);
%! assert (size (m.CoefficientCovariance), [2, 2]);
%! assert (m.CoefficientNames{1}, 'x1');
%! assert (m.CoefficientNames{2}, 'x2');
%! assert (! any (strcmp (m.CoefficientNames, '(Intercept)')));
%! assert (height (m.Diagnostics), 20);

%!test
%! ## removing both predictors leaves only the intercept
%! m = removeTerms (mdl, 'x1 + x2');
%! assert (m.NumCoefficients, 1);
%! assert (m.NumEstimatedCoefficients, 1);
%! assert (m.DFE, 19);
%! assert (m.NumObservations, 20);
%! assert (m.Formula.HasIntercept, true);
%! assert (m.Formula.LinearPredictor, '1');
%! assert (m.Coefficients.Estimate(1), -5.5900177811558, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.SSE, 583.910420002346, -1e-8);
%! assert (m.SST, 583.910420002346, -1e-8);
%! assert (m.SSR, 0, 1e-20);
%! assert (size (m.CoefficientCovariance), [1, 1]);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (height (m.Diagnostics), 20);
%! assert (sum (m.Diagnostics.Leverage), 1, 1e-10);
%! assert (m.Residuals.Raw, y - m.Fitted, 1e-10);

%!test
%! warning ('off', 'all');
%! lastwarn ('');
%! m = removeTerms (mdl, 'x1:x2');
%! wmsg = lastwarn ();
%! warning ('on', 'all');
%! assert (! isempty (strfind (wmsg, 'No specified terms appear in the model')));
%! assert (m.NumCoefficients, mdl.NumCoefficients);
%! assert (m.NumEstimatedCoefficients, mdl.NumEstimatedCoefficients);
%! assert (m.DFE, mdl.DFE);
%! assert (m.SSE, mdl.SSE, 1e-15);
%! assert (m.SSR, mdl.SSR, 1e-15);
%! assert (m.SST, mdl.SST, 1e-15);
%! assert (m.RMSE, mdl.RMSE, 1e-15);
%! assert (m.Coefficients.Estimate, mdl.Coefficients.Estimate, 1e-15);
%! assert (m.Coefficients.SE, mdl.Coefficients.SE, 1e-15);
%! assert (m.CoefficientCovariance, mdl.CoefficientCovariance, 1e-15);
%! assert (isequal (m.CoefficientNames, mdl.CoefficientNames));
%! assert (m.Formula.LinearPredictor, mdl.Formula.LinearPredictor);
%! assert (m.Formula.HasIntercept, mdl.Formula.HasIntercept);

%!test
%! m = removeTerms (mdl, [0 1 0]);
%! assert (m.NumCoefficients, 2);
%! assert (m.NumEstimatedCoefficients, 2);
%! assert (m.DFE, 18);
%! assert (m.NumObservations, 20);
%! assert (m.Coefficients.Estimate(1), 3.88470469761717, -1e-8);
%! assert (m.Coefficients.Estimate(2), -18.047090435758, -1e-8);
%! assert (m.Coefficients.SE(2), 1.19086428900602, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.Coefficients.tStat(2), -15.1546155194741, -1e-8);
%! assert (m.SSE, 42.4383708132815, -1e-8);
%! assert (m.SST, 583.910420002346, -1e-8);
%! assert (m.Rsquared.Ordinary, 0.927320408474452, -1e-8);
%! assert (size (m.CoefficientCovariance), [2, 2]);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.Formula.HasIntercept, true);
%! assert (m.Formula.LinearPredictor, '1 + x1');
%! assert (height (m.Diagnostics), 20);
%! assert (sum (m.Diagnostics.Leverage), 2, 1e-10);

%!test
%! ## auto-padded matrix [0 1] gives identical result to [0 1 0]
%! m = removeTerms (mdl, [0 1]);
%! assert (m.NumCoefficients, 2);
%! assert (m.NumEstimatedCoefficients, 2);
%! assert (m.DFE, 18);
%! assert (m.Coefficients.Estimate(1), 3.88470469761717, -1e-8);
%! assert (m.Coefficients.Estimate(2), -18.047090435758, -1e-8);
%! assert (m.Coefficients.SE(2), 1.19086428900602, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.SSE, 42.4383708132815, -1e-8);
%! assert (m.SST, 583.910420002346, -1e-8);
%! assert (m.Rsquared.Ordinary, 0.927320408474452, -1e-8);
%! assert (size (m.CoefficientCovariance), [2, 2]);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.Formula.LinearPredictor, '1 + x1');
%! assert (height (m.Diagnostics), 20);
%! assert (sum (m.Diagnostics.Leverage), 2, 1e-10);

%!test
%! ## multi-row matrix removes two terms same as the string form
%! Xh = [7 26 6 60; 1 29 15 52; 11 56 8 20; 11 31 8 47; 7 52 6 33; ...
%!        11 55 9 22; 3 71 17 6; 1 31 22 44; 2 54 18 22; 21 47 4 26; ...
%!        1 40 23 34; 11 66 9 12; 10 68 8 12];
%! yh = [78.5;74.3;104.3;87.6;95.9;109.2;102.7;72.5;93.1;115.9;83.8;113.3;109.4];
%! m = removeTerms (fitlm (Xh, yh), [0 0 1 0 0; 0 0 0 1 0]);
%! assert (m.NumCoefficients, 3);
%! assert (m.NumEstimatedCoefficients, 3);
%! assert (m.DFE, 10);
%! assert (m.NumObservations, 13);
%! assert (m.Coefficients.Estimate(1), 52.577348882089481, -1e-8);
%! assert (m.Coefficients.Estimate(2), 1.468305742215555, -1e-8);
%! assert (m.Coefficients.Estimate(3), 0.662250491274645, -1e-8);
%! assert (m.Coefficients.SE(1), 2.286174334503340, -1e-8);
%! assert (m.Coefficients.SE(2), 0.121300923606266, -1e-8);
%! assert (m.Coefficients.SE(3), 0.045854721468522, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.Coefficients.tStat(2), 12.104654264476748, -1e-8);
%! assert (m.Coefficients.tStat(3), 14.442362096327519, -1e-8);
%! assert (m.Coefficients.pValue(2), 2.692212179685427e-07, -1e-8);
%! assert (m.SSE, 57.904483176113658, -1e-8);
%! assert (m.Rsquared.Ordinary, 0.978678374535632, -1e-8);
%! assert (m.Rsquared.Adjusted, 0.974414049442758, -1e-8);
%! assert (size (m.CoefficientCovariance), [3, 3]);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.CoefficientNames{3}, 'x2');
%! assert (m.Formula.HasIntercept, true);
%! assert (m.Formula.LinearPredictor, '1 + x1 + x2');
%! assert (height (m.Diagnostics), 13);
%! assert (sum (m.Diagnostics.Leverage), 3, 1e-10);
%! assert (m.Residuals.Raw, yh - m.Fitted, 1e-10);

%!test
%! ## observation weights carry through to the refitted model
%! w = (1:n)' / sum (1:n);
%! mw = fitlm (X, y, 'Weights', w);
%! m = removeTerms (mw, 'x2');
%! assert (m.NumCoefficients, 2);
%! assert (m.NumEstimatedCoefficients, 2);
%! assert (m.DFE, 18);
%! assert (m.NumObservations, 20);
%! assert (m.Coefficients.Estimate(1), 6.29263960898714, -1e-8);
%! assert (m.Coefficients.Estimate(2), -21.5708976231287, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.SSE, 1.41763159723151, -1e-8);
%! assert (size (m.CoefficientCovariance), [2, 2]);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.Formula.HasIntercept, true);
%! assert (m.Formula.LinearPredictor, '1 + x1');
%! assert (m.ObservationInfo.Weights, w, 1e-15);
%! assert (sum (m.ObservationInfo.Weights), 1, 1e-12);
%! assert (height (m.Diagnostics), 20);
%! assert (sum (m.Diagnostics.Leverage), 2, 1e-10);
%! assert (m.SSE != removeTerms (mdl, 'x2').SSE);

%!test
%! ## excluded rows are preserved and reduce effective sample size
%! me = fitlm (X, y, 'Exclude', [1, 3]);
%! m = removeTerms (me, 'x2');
%! assert (m.NumObservations, 18);
%! assert (m.DFE, 16);
%! assert (m.NumCoefficients, 2);
%! assert (m.NumEstimatedCoefficients, 2);
%! assert (m.Coefficients.Estimate(1), 4.96609542902066, -1e-8);
%! assert (m.Coefficients.Estimate(2), -19.5618050042778, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (size (m.CoefficientCovariance), [2, 2]);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.Formula.LinearPredictor, '1 + x1');
%! assert (m.ObservationInfo.Excluded(1), true);
%! assert (m.ObservationInfo.Excluded(3), true);
%! assert (m.ObservationInfo.Excluded(2), false);
%! assert (m.ObservationInfo.Missing(1), false);
%! assert (height (m.Diagnostics), 20);
%! assert (isnan (m.Fitted(1)));
%! assert (isnan (m.Fitted(3)));
%! assert (isfinite (m.Fitted(2)));

%!test
%! ## removing x2 from a no-intercept model gives one slope term
%! mni = fitlm (X, y, 'Intercept', false);
%! m = removeTerms (mni, 'x2');
%! assert (m.NumCoefficients, 1);
%! assert (m.NumEstimatedCoefficients, 1);
%! assert (m.DFE, 19);
%! assert (m.NumObservations, 20);
%! assert (m.Formula.HasIntercept, false);
%! assert (m.Formula.LinearPredictor, 'x1');
%! assert (m.Coefficients.Estimate(1), -12.362156731928, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.SSE, 112.371951585499, -1e-8);
%! assert (size (m.CoefficientCovariance), [1, 1]);
%! assert (m.CoefficientCovariance(1,1) > 0);
%! assert (m.CoefficientNames{1}, 'x1');
%! assert (! any (strcmp (m.CoefficientNames, '(Intercept)')));
%! assert (height (m.Diagnostics), 20);
%! assert (all (isfinite (m.Fitted)));
%! assert (sum (m.Diagnostics.Leverage), 1, 1e-10);

%!test
%! ## removing the interaction term recovers the plain linear model
%! mi = fitlm (X, y, 'interactions');
%! m = removeTerms (mi, 'x1:x2');
%! assert (m.NumCoefficients, 3);
%! assert (m.NumEstimatedCoefficients, 3);
%! assert (m.DFE, 17);
%! assert (m.NumObservations, 20);
%! assert (m.Coefficients.Estimate(1), 0.116188677790207, -1e-8);
%! assert (m.Coefficients.Estimate(2), 2.508451490570863, -1e-8);
%! assert (m.Coefficients.Estimate(3), -0.978835329825186, -1e-8);
%! assert (m.Coefficients.SE(1), 0.112185831, -1e-7);
%! assert (m.Coefficients.SE(2), 0.4920818186, -1e-7);
%! assert (m.Coefficients.SE(3), 0.02276108523, -1e-7);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.Coefficients.tStat(1), 1.035680502, -1e-6);
%! assert (m.Coefficients.tStat(2), 5.097630913, -1e-6);
%! assert (m.Coefficients.tStat(3), -43.00477415, -1e-6);
%! assert (m.SSE, 0.386545331386823, -1e-8);
%! assert (size (m.CoefficientCovariance), [3, 3]);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.CoefficientNames{3}, 'x2');
%! assert (m.Formula.HasIntercept, true);
%! assert (m.Formula.LinearPredictor, '1 + x1 + x2');
%! assert (height (m.Diagnostics), 20);
%! assert (sum (m.Diagnostics.Leverage), 3, 1e-10);

%!test
%! ## removing a quadratic term refits on the remaining terms
%! mq = fitlm (X, y, 'quadratic');
%! m = removeTerms (mq, 'x2^2');
%! assert (m.NumEstimatedCoefficients, 4);
%! assert (m.DFE, 16);
%! assert (m.SSE, 0.383859187927621, -1e-8);
%! assert (size (m.CoefficientCovariance, 1), m.NumCoefficients);
%! assert (size (m.CoefficientCovariance, 2), m.NumCoefficients);
%! assert (m.Formula.HasIntercept, true);
%! assert (height (m.Diagnostics), 20);
%! assert (sum (m.Diagnostics.Leverage), 4, 1e-10);
%! assert (m.SSE >= mq.SSE);
%! assert (! any (strcmp (m.CoefficientNames, 'x2^2')));

%!test
%! ## star notation removes main effects and interaction in one call
%! mi = fitlm (X, y, 'interactions');
%! m = removeTerms (mi, 'x1*x2');
%! assert (m.NumCoefficients, 1);
%! assert (m.NumEstimatedCoefficients, 1);
%! assert (m.DFE, 19);
%! assert (m.NumObservations, 20);
%! assert (m.Formula.HasIntercept, true);
%! assert (m.Formula.LinearPredictor, '1');
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.Coefficients.Estimate(1), -5.5900177811558, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.SSE, 583.910420002346, -1e-8);
%! assert (m.SST, 583.910420002346, -1e-8);
%! assert (m.SSR, 0, 1e-20);
%! assert (size (m.CoefficientCovariance), [1, 1]);
%! assert (height (m.Diagnostics), 20);
%! assert (sum (m.Diagnostics.Leverage), 1, 1e-10);

%!test
%! ## 3-predictor model: removing one term matches a direct two-predictor fit
%! X3 = [X, sin((1:n)' * pi / n)];
%! y3 = X3 * [3; -1; 2] + 0.1 * cos ((1:n)' * pi / 7);
%! m = removeTerms (fitlm (X3, y3), 'x3');
%! r = fitlm (X, y3);
%! assert (m.NumCoefficients, 3);
%! assert (m.NumEstimatedCoefficients, 3);
%! assert (m.DFE, 17);
%! assert (m.NumObservations, 20);
%! assert (m.Coefficients.Estimate, r.Coefficients.Estimate, 1e-10);
%! assert (m.Coefficients.SE, r.Coefficients.SE, 1e-10);
%! assert (m.Coefficients.tStat, r.Coefficients.tStat, 1e-10);
%! assert (m.Coefficients.pValue, r.Coefficients.pValue, 1e-10);
%! assert (m.SSE, r.SSE, 1e-12);
%! assert (m.SSR, r.SSR, 1e-12);
%! assert (m.SST, r.SST, 1e-12);
%! assert (m.MSE, r.MSE, 1e-12);
%! assert (m.RMSE, r.RMSE, 1e-12);
%! assert (m.Rsquared.Ordinary, r.Rsquared.Ordinary, 1e-12);
%! assert (m.Rsquared.Adjusted, r.Rsquared.Adjusted, 1e-12);
%! assert (m.CoefficientCovariance, r.CoefficientCovariance, 1e-12);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.CoefficientNames{3}, 'x2');
%! assert (m.Formula.HasIntercept, true);
%! assert (m.Formula.LinearPredictor, '1 + x1 + x2');
%! assert (height (m.Diagnostics), 20);
%! assert (sum (m.Diagnostics.Leverage), 3, 1e-10);

%!test
%! ## 3-predictor model: removing two terms matches a direct one-predictor fit
%! X3 = [X, sin((1:n)' * pi / n)];
%! y3 = X3 * [3; -1; 2] + 0.1 * cos ((1:n)' * pi / 7);
%! m = removeTerms (fitlm (X3, y3), 'x2 + x3');
%! r = fitlm (X(:,1), y3);
%! assert (m.NumCoefficients, 2);
%! assert (m.NumEstimatedCoefficients, 2);
%! assert (m.DFE, 18);
%! assert (m.NumObservations, 20);
%! assert (m.Coefficients.Estimate, r.Coefficients.Estimate, 1e-10);
%! assert (m.Coefficients.SE, r.Coefficients.SE, 1e-10);
%! assert (m.Coefficients.tStat, r.Coefficients.tStat, 1e-10);
%! assert (m.Coefficients.pValue, r.Coefficients.pValue, 1e-10);
%! assert (m.SSE, r.SSE, 1e-12);
%! assert (m.SST, r.SST, 1e-12);
%! assert (m.MSE, r.MSE, 1e-12);
%! assert (m.RMSE, r.RMSE, 1e-12);
%! assert (m.Rsquared.Ordinary, r.Rsquared.Ordinary, 1e-12);
%! assert (m.Rsquared.Adjusted, r.Rsquared.Adjusted, 1e-12);
%! assert (m.CoefficientCovariance, r.CoefficientCovariance, 1e-12);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.Formula.HasIntercept, true);
%! assert (m.Formula.LinearPredictor, '1 + x1');
%! assert (height (m.Diagnostics), 20);
%! assert (sum (m.Diagnostics.Leverage), 2, 1e-10);

%!test
%! Xh = [7 26 6 60; 1 29 15 52; 11 56 8 20; 11 31 8 47; 7 52 6 33; ...
%!        11 55 9 22; 3 71 17 6; 1 31 22 44; 2 54 18 22; 21 47 4 26; ...
%!        1 40 23 34; 11 66 9 12; 10 68 8 12];
%! yh = [78.5;74.3;104.3;87.6;95.9;109.2;102.7;72.5;93.1;115.9;83.8;113.3;109.4];
%! m = removeTerms (fitlm (Xh, yh), 'x3 + x4');
%! r = removeTerms (removeTerms (fitlm (Xh, yh), 'x4'), 'x3');
%! assert (r.NumCoefficients, 3);
%! assert (r.DFE, 10);
%! assert (r.SSE, m.SSE, 1e-12);
%! assert (r.SSR, m.SSR, 1e-12);
%! assert (r.SST, m.SST, 1e-12);
%! assert (r.RMSE, m.RMSE, 1e-12);
%! assert (r.Coefficients.Estimate, m.Coefficients.Estimate, 1e-10);
%! assert (r.Coefficients.SE, m.Coefficients.SE, 1e-10);
%! assert (r.Coefficients.tStat, m.Coefficients.tStat, 1e-10);
%! assert (r.Coefficients.pValue, m.Coefficients.pValue, 1e-10);
%! assert (r.CoefficientCovariance, m.CoefficientCovariance, 1e-12);
%! assert (isequal (r.CoefficientNames, m.CoefficientNames));
%! assert (r.Rsquared.Ordinary, m.Rsquared.Ordinary, 1e-12);
%! assert (r.Rsquared.Adjusted, m.Rsquared.Adjusted, 1e-12);
%! assert (r.Formula.LinearPredictor, m.Formula.LinearPredictor);

%!test
%! m = removeTerms (mdl, [0 0 0]);
%! r = removeTerms (mdl, '1');
%! assert (m.NumCoefficients, 2);
%! assert (m.Formula.HasIntercept, false);
%! assert (m.Formula.LinearPredictor, 'x1 + x2');
%! assert (m.Coefficients.Estimate(1), 2.96142161317611, -1e-8);
%! assert (m.Coefficients.Estimate(2), -0.997248749443286, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.SSE, r.SSE, 1e-15);
%! assert (m.SSR, r.SSR, 1e-15);
%! assert (m.SST, r.SST, 1e-15);
%! assert (m.Coefficients.Estimate, r.Coefficients.Estimate, 1e-15);
%! assert (m.Coefficients.SE, r.Coefficients.SE, 1e-15);
%! assert (m.CoefficientCovariance, r.CoefficientCovariance, 1e-15);
%! assert (isequal (m.CoefficientNames, r.CoefficientNames));
%! assert (height (m.Diagnostics), 20);
%! assert (sum (m.Diagnostics.Leverage), 2, 1e-10);

%!test
%! ## removing a categorical predictor drops all its indicator variables at once
%! Xc = [1;1;1;2;2;2;3;3;3];
%! yc = [2.1;2.3;1.9; 4.1;3.9;4.2; 6.3;5.8;6.1];
%! mc = fitlm (Xc, yc, 'linear', 'CategoricalVars', 1);
%! m = removeTerms (mc, 'x1');
%! assert (m.NumCoefficients, 1);
%! assert (m.NumEstimatedCoefficients, 1);
%! assert (m.DFE, 8);
%! assert (m.NumObservations, 9);
%! assert (m.Formula.HasIntercept, true);
%! assert (m.Formula.LinearPredictor, '1');
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.Coefficients.Estimate(1), 4.07777777777778, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.SSR, 0, 1e-20);
%! assert (size (m.CoefficientCovariance), [1, 1]);
%! assert (height (m.Diagnostics), 9);
%! assert (all (isfinite (m.Fitted)));
%! assert (sum (m.Diagnostics.Leverage), 1, 1e-10);
%! assert (m.Residuals.Raw, yc - m.Fitted, 1e-10);

%!test
%! ## matrix row removes x4 from hald leaving intercept plus x1 x2 x3
%! Xh = [7 26 6 60; 1 29 15 52; 11 56 8 20; 11 31 8 47; 7 52 6 33; ...
%!        11 55 9 22; 3 71 17 6; 1 31 22 44; 2 54 18 22; 21 47 4 26; ...
%!        1 40 23 34; 11 66 9 12; 10 68 8 12];
%! yh = [78.5;74.3;104.3;87.6;95.9;109.2;102.7;72.5;93.1;115.9;83.8;113.3;109.4];
%! m = removeTerms (fitlm (Xh, yh), [0 0 0 1 0]);
%! assert (m.NumCoefficients, 4);
%! assert (m.NumEstimatedCoefficients, 4);
%! assert (m.DFE, 9);
%! assert (m.NumObservations, 13);
%! assert (m.Coefficients.Estimate(1), 48.1936343180437, -1e-8);
%! assert (m.Coefficients.Estimate(2), 1.69589016748479, -1e-8);
%! assert (m.Coefficients.Estimate(3), 0.656914878270554, -1e-8);
%! assert (m.Coefficients.Estimate(4), 0.250017606680009, -1e-8);
%! assert (m.Coefficients.tStat, m.Coefficients.Estimate ./ m.Coefficients.SE, 1e-10);
%! assert (m.SSE, 48.1106140726532, -1e-8);
%! assert (size (m.CoefficientCovariance), [4, 4]);
%! assert (m.CoefficientNames{1}, '(Intercept)');
%! assert (m.CoefficientNames{2}, 'x1');
%! assert (m.CoefficientNames{3}, 'x2');
%! assert (m.CoefficientNames{4}, 'x3');
%! assert (m.Formula.HasIntercept, true);
%! assert (m.Formula.LinearPredictor, '1 + x1 + x2 + x3');
%! assert (height (m.Diagnostics), 13);
%! assert (all (isfinite (m.Fitted)));
%! assert (sum (m.Diagnostics.Leverage), 4, 1e-10);
%! assert (m.Residuals.Raw, yh - m.Fitted, 1e-10);

%!test
%! ## default call creates a histogram with correct bin count and density
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl);
%! xd = get (h(1), 'XData');
%! yd = get (h(1), 'YData');
%! r = mdl.Residuals.Raw(! isnan (mdl.Residuals.Raw));
%! bw = xd(3,1) - xd(1,1);
%! assert (numel (h), 1);
%! assert (get (h(1), 'type'), 'patch');
%! assert (size (xd, 2) > 0);
%! assert (sum (yd(2,:)) * bw, 1, 1e-10);
%! assert (all (yd(1,:) == 0) && all (yd(4,:) == 0));
%! assert (get (get (ax, 'xlabel'), 'string'), 'Residuals');
%! assert (get (get (ax, 'ylabel'), 'string'), 'Probability density');
%! assert (get (get (ax, 'title'), 'string'), 'Histogram of residuals');
%! close (fig);

%!test
%! ## histogram bar color changes when FaceColor is passed
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'histogram', 'FaceColor', [0 1 0]);
%! assert (get (h(1), 'FaceColor'), [0 1 0], 1e-10);
%! close (fig);

%!test
%! ## fitted plot shows residuals against fitted values with a zero reference line
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'fitted');
%! assert (numel (h), 2);
%! assert (get (h(1), 'XData'), mdl.Fitted', 1e-15);
%! assert (get (h(1), 'YData'), mdl.Residuals.Raw', 1e-15);
%! assert (get (h(1), 'LineStyle'), 'none');
%! assert (get (h(1), 'Marker'), 'x');
%! assert (get (h(2), 'YData'), [0 0]);
%! assert (get (h(2), 'LineStyle'), ':');
%! assert (get (get (ax, 'xlabel'), 'string'), 'Fitted values');
%! assert (get (get (ax, 'ylabel'), 'string'), 'Residuals');
%! assert (get (get (ax, 'title'), 'string'), 'Plot of residuals vs. fitted values');
%! close (fig);

%!test
%! ## custom color applies to data points but leaves the reference line unchanged
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'fitted', 'Color', [1 0 0]);
%! assert (get (h(1), 'Color'), [1 0 0], 1e-10);
%! assert (get (h(2), 'Color'), [0.8510 0.8510 0.8510], 1e-4);
%! assert (get (h(2), 'LineStyle'), ':');
%! close (fig);

%!test
%! ## excluded rows appear as gaps in the fitted plot
%! me = fitlm (X, y, 'Exclude', [3, 8]);
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, me, 'fitted');
%! yd = get (h(1), 'YData');
%! assert (numel (yd), 20);
%! assert (isnan (yd(3)));
%! assert (isnan (yd(8)));
%! assert (! isnan (yd(1)));
%! close (fig);

%!test
%! ## case order plot covers all rows and shows gaps where rows were excluded
%! me = fitlm (X, y, 'Exclude', [2, 5]);
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, me, 'caseorder');
%! xd = get (h(1), 'XData');
%! yd = get (h(1), 'YData');
%! assert (xd, 1:20);
%! assert (isnan (yd(2)));
%! assert (isnan (yd(5)));
%! assert (! isnan (yd(1)));
%! assert (get (h(2), 'YData'), [0 0]);
%! assert (get (get (ax, 'xlabel'), 'string'), 'Row number');
%! assert (get (get (ax, 'ylabel'), 'string'), 'Residuals');
%! assert (get (get (ax, 'title'), 'string'), 'Case order plot of residuals');
%! close (fig);

%!test
%! ## lagged plot shows each residual against the previous one with two reference lines
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'lagged');
%! r = mdl.Residuals.Raw;
%! assert (numel (h), 3);
%! assert (get (h(1), 'XData'), r(1:end-1)', 1e-15);
%! assert (get (h(1), 'YData'), r(2:end)', 1e-15);
%! assert (get (h(2), 'YData'), [0 0]);
%! assert (get (h(3), 'XData'), [0 0]);
%! assert (get (h(2), 'LineStyle'), ':');
%! assert (get (get (ax, 'xlabel'), 'string'), 'Residual(t-1)');
%! assert (get (get (ax, 'ylabel'), 'string'), 'Residual(t)');
%! assert (get (get (ax, 'title'), 'string'), 'Plot of residuals vs. lagged residuals');
%! close (fig);

%!test
%! ## probability plot uses sorted active residuals as its data
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'probability');
%! r_s = sort (mdl.Residuals.Raw(! isnan (mdl.Residuals.Raw)));
%! assert (get (h(1), 'XData'), r_s', 1e-15);
%! assert (get (get (ax, 'xlabel'), 'string'), 'Residuals');
%! assert (get (get (ax, 'ylabel'), 'string'), 'Probability');
%! assert (get (get (ax, 'title'), 'string'), 'Normal probability plot of residuals');
%! close (fig);

%!test
%! ## observed plot connects each point to the reference line with a vertical segment
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'observed');
%! obs = mdl.Variables{:, mdl.ResponseName};
%! assert (numel (h), 3);
%! assert (get (h(1), 'XData'), mdl.Fitted', 1e-15);
%! assert (get (h(1), 'YData'), obs', 1e-15);
%! assert (isequal (get (h(2), 'XData'), get (h(2), 'YData')));
%! xd3 = get (h(3), 'XData');
%! assert (numel (xd3), 3 * mdl.NumObservations);
%! assert (sum (isnan (xd3)), mdl.NumObservations);
%! assert (get (get (ax, 'xlabel'), 'string'), 'Fitted values');
%! assert (get (get (ax, 'ylabel'), 'string'), 'Observed response values');
%! assert (get (get (ax, 'title'), 'string'), 'Plot of observed vs. fitted values');
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
%! assert (numel (h), 2);
%! assert (get (h(1), 'XData'), x_sym', 1e-15);
%! assert (get (h(1), 'YData'), y_sym', 1e-15);
%! assert (isequal (get (h(2), 'XData'), get (h(2), 'YData')));
%! assert (get (h(2), 'LineStyle'), ':');
%! assert (get (get (ax, 'xlabel'), 'string'), 'Lower tail');
%! assert (get (get (ax, 'ylabel'), 'string'), 'Upper tail');
%! close (fig);

%!test
%! ## switching to pearson residuals changes plotted values but not x positions
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'fitted', 'ResidualType', 'pearson');
%! assert (get (h(1), 'YData'), mdl.Residuals.Pearson', 1e-15);
%! assert (get (h(1), 'XData'), mdl.Fitted', 1e-15);
%! assert (! isequal (get (h(1), 'YData'), mdl.Residuals.Raw'));
%! close (fig);

%!test
%! ## standardized and studentized residuals produce different values
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h1 = plotResiduals (ax, mdl, 'caseorder', 'ResidualType', 'standardized');
%! h2 = plotResiduals (ax, mdl, 'caseorder', 'ResidualType', 'studentized');
%! assert (get (h1(1), 'YData'), mdl.Residuals.Standardized', 1e-15);
%! assert (get (h2(1), 'YData'), mdl.Residuals.Studentized', 1e-15);
%! assert (! isequal (get (h1(1), 'YData'), get (h2(1), 'YData')));
%! close (fig);

%!test
%! ## marker style and size apply to data points but not to the reference line
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mdl, 'fitted', 'Marker', 's', 'MarkerSize', 10);
%! assert (get (h(1), 'Marker'), 's');
%! assert (get (h(1), 'MarkerSize'), 10);
%! assert (get (h(2), 'Marker'), 'none');
%! close (fig);

%!test
%! ## weighted model residuals differ from unweighted residuals in the fitted plot
%! mw = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotResiduals (ax, mw, 'fitted');
%! assert (get (h(1), 'YData'), mw.Residuals.Raw', 1e-15);
%! assert (! isequal (get (h(1), 'YData'), mdl.Residuals.Raw'));
%! close (fig);

%!test
%! ## calling without an axes handle plots into the current axes
%! fig = figure ('visible', 'off');
%! h = plotResiduals (mdl, 'fitted');
%! assert (isgraphics (get (h(1), 'Parent'), 'axes'));
%! assert (isequal (get (h(1), 'Parent'), gca ()));
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl);
%! yd = get (h(1), 'YData');
%! assert (numel (h), 2);
%! assert (get (h(1), 'XData'), 1:n);
%! assert (yd(1), 0.370779220779221, -1e-10);
%! assert (yd(2), 0.245283663704716, -1e-10);
%! assert (yd(3), 0.164718614718615, -1e-10);
%! assert (yd(4), 0.118147641831852, -1e-10);
%! assert (yd(5), 0.0960013670539986, -1e-10);
%! assert (yd(6), 0.0900774663932558, -1e-10);
%! assert (yd(7), 0.0935406698564593, -1e-10);
%! assert (yd(8), 0.100922761449077, -1e-10);
%! assert (yd(9), 0.108122579175211, -1e-10);
%! assert (yd(10), 0.112406015037594, -1e-10);
%! assert (yd(11), 0.112406015037594, -1e-10);
%! assert (yd(12), 0.108122579175211, -1e-10);
%! assert (yd(13), 0.100922761449077, -1e-10);
%! assert (yd(14), 0.0935406698564592, -1e-10);
%! assert (yd(15), 0.0900774663932559, -1e-10);
%! assert (yd(16), 0.0960013670539986, -1e-10);
%! assert (yd(17), 0.118147641831852, -1e-10);
%! assert (yd(18), 0.164718614718615, -1e-10);
%! assert (yd(19), 0.245283663704716, -1e-10);
%! assert (yd(20), 0.370779220779221, -1e-10);
%! assert (get (h(2), 'YData'), [0.3, 0.3], 1e-12);
%! assert (get (h(2), 'XData'), [0, n]);
%! assert (get (h(2), 'LineStyle'), ':');
%! assert (get (get (ax, 'xlabel'), 'string'), 'Row number');
%! assert (get (get (ax, 'ylabel'), 'string'), 'Leverage');
%! assert (get (get (ax, 'title'), 'string'), 'Case order plot of leverage');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'leverage', 'Color', [1 0 0]);
%! assert (get (h(1), 'Color'), [1 0 0], 1e-10);
%! assert (get (h(2), 'Color'), [0.8510 0.8510 0.8660], 1e-4);
%! assert (get (h(2), 'LineStyle'), ':');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'cookd');
%! yd = get (h(1), 'YData');
%! assert (numel (h), 2);
%! assert (yd(1), 0.078517048682575, -1e-8);
%! assert (yd(2), 0.077211407930332, -1e-8);
%! assert (yd(3), 0.001953301452841, -1e-7);
%! assert (get (h(2), 'YData'), [0.1668641787, 0.1668641787], -1e-8);
%! assert (get (h(2), 'XData'), [0, n]);
%! assert (get (get (ax, 'ylabel'), 'string'), 'Cook''s distance');
%! assert (get (get (ax, 'title'), 'string'), 'Case order plot of Cook''s distance');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'covratio');
%! yd = get (h(1), 'YData');
%! yv = get (h(2), 'YData');
%! xv = get (h(2), 'XData');
%! assert (numel (h), 2);
%! assert (yd(1), 1.774933177, -1e-8);
%! assert (yd(2), 1.397661919, -1e-8);
%! assert (yd(3), 1.428481535, -1e-8);
%! assert (numel (xv), 5);
%! assert (sum (isnan (xv)), 1);
%! assert (yv(1), 0.55, 1e-12);
%! assert (yv(2), 0.55, 1e-12);
%! assert (yv(4), 1.45, 1e-12);
%! assert (yv(5), 1.45, 1e-12);
%! assert (get (get (ax, 'ylabel'), 'string'), 'Covariance ratio');
%! assert (get (get (ax, 'title'), 'string'), 'Case order plot of covariance ratio');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'dfbetas');
%! p = mdl.NumCoefficients;
%! yv = get (h(p+1), 'YData');
%! xv = get (h(p+1), 'XData');
%! assert (numel (h), p + 1);
%! assert (numel (get (h(1), 'YData')), n);
%! assert (numel (get (h(2), 'YData')), n);
%! assert (numel (get (h(3), 'YData')), n);
%! assert (numel (xv), 5);
%! assert (sum (isnan (xv)), 1);
%! assert (yv(1), -0.6708203932, -1e-8);
%! assert (yv(end), 0.6708203932, -1e-8);
%! assert (get (get (ax, 'ylabel'), 'string'), 'Scaled change in coefficients');
%! assert (get (get (ax, 'title'), 'string'), 'Case order plot of scaled change in coefficients');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'dfbetas', 'Color', [1 0 0]);
%! p = mdl.NumCoefficients;
%! for k = 1:p
%! assert (get (h(k), 'Color'), [1 0 0], 1e-10);
%! endfor
%! assert (get (h(p+1), 'Color'), [0.8510 0.8510 0.8660], 1e-4);
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'dffits');
%! yd = get (h(1), 'YData');
%! yv = get (h(2), 'YData');
%! xv = get (h(2), 'XData');
%! assert (numel (h), 2);
%! assert (yd(1), 0.476480465355394, -1e-8);
%! assert (yd(2), 0.477020506700835, -1e-8);
%! assert (yd(3), -0.074329411030064, -1e-7);
%! assert (sum (isnan (xv)), 1);
%! assert (yv(1), -0.7745966692, -1e-8);
%! assert (yv(end), 0.7745966692, -1e-8);
%! assert (get (get (ax, 'ylabel'), 'string'), 'Scaled change in fit');
%! assert (get (get (ax, 'title'), 'string'), 'Case order plot of scaled change in fit');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 's2_i');
%! yd = get (h(1), 'YData');
%! assert (numel (h), 2);
%! assert (yd(1), 0.02359100986, -1e-8);
%! assert (yd(2), 0.02314622330, -1e-8);
%! assert (yd(3), 0.02411685408, -1e-8);
%! assert (get (h(2), 'YData'), [0.02273796067, 0.02273796067], -1e-8);
%! assert (get (h(2), 'XData'), [0, n]);
%! assert (get (get (ax, 'ylabel'), 'string'), 'Leave-one-out variance');
%! assert (get (get (ax, 'title'), 'string'), 'Case order plot of leave-one-out variance');
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotDiagnostics (ax, mdl, 'contour');
%! yd = get (h(1), 'YData');
%! xd = get (h(1), 'XData');
%! assert (numel (h), 2);
%! assert (xd(1), 0.3707792208, -1e-8);
%! assert (xd(2), 0.2452836637, -1e-8);
%! assert (yd(1), 0.07562471113, -1e-8);
%! assert (yd(2), 0.11059272450, -1e-8);
%! assert (get (h(1), 'LineStyle'), 'none');
%! assert (get (h(1), 'Marker'), 'x');
%! assert (isgraphics (h(2)));
%! assert (get (get (ax, 'xlabel'), 'string'), 'Leverage');
%! assert (get (get (ax, 'ylabel'), 'string'), 'Residual');
%! assert (get (get (ax, 'title'), 'string'), 'Cook''s distance factorization');
%! close (fig);

%!test
%! me = fitlm (X, y, 'Exclude', [2, 7]);
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h_ex = plotDiagnostics (ax, me, 'cookd');
%! h_un = plotDiagnostics (ax, mdl, 'cookd');
%! ref_ex = get (h_ex(2), 'YData');
%! ref_un = get (h_un(2), 'YData');
%! assert (! isequal (ref_ex, ref_un));
%! assert (ref_ex(1), 3 * mean (me.Diagnostics.CooksDistance, 'omitnan'), 1e-12);
%! close (fig);

%!test
%! mw = fitlm (X, y, 'Weights', (1:n)' / sum (1:n));
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! hw = plotDiagnostics (ax, mw, 'leverage');
%! hu = plotDiagnostics (ax, mdl, 'leverage');
%! ydw = get (hw(1), 'YData');
%! ydu = get (hu(1), 'YData');
%! assert (ydw(1) != ydu(1));
%! assert (! isequal (ydw, ydu));
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! h = plotDiagnostics (mdl);
%! assert (isequal (get (h(1), 'Parent'), gca ()));
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
%! assert (numel (h), 3);
%! assert (xd1(1), 2.38302891604232, -1e-10);
%! assert (xd1(2), -19.5277648300125, -1e-10);
%! assert (yd1, [1 2]);
%! assert (xd2(1), 1.39673712385796, -1e-10);
%! assert (xd2(2), 3.36932070822668, -1e-10);
%! assert (yd2, [1 1]);
%! assert (xd3(1), -20.4857975891918, -1e-10);
%! assert (xd3(2), -18.5697320708331, -1e-10);
%! assert (yd3, [2 2]);
%! assert (get (h(1), 'Color'), [0.1490 0.5490 0.8660], 1e-4);
%! assert (get (h(2), 'Color'), [0.1490 0.5490 0.8660], 1e-4);
%! assert (get (h(3), 'Color'), [0.1490 0.5490 0.8660], 1e-4);
%! assert (get (h(1), 'Marker'), 'o');
%! assert (get (h(1), 'LineStyle'), 'none');
%! assert (get (h(2), 'LineStyle'), '-');
%! assert (get (h(2), 'Marker'), 'none');
%! assert (get (h(3), 'LineStyle'), '-');
%! assert (get (h(3), 'Marker'), 'none');
%! assert (mean (xd2), xd1(1), 1e-10);
%! assert (mean (xd3), xd1(2), 1e-10);
%! assert (get (get (ax, 'xlabel'), 'string'), 'Main Effect');
%! assert (get (get (ax, 'ylabel'), 'string'), '');
%! assert (get (get (ax, 'title'), 'string'), 'Main Effects Plot');
%! assert (get (ax, 'YTick'), [1 2]);
%! assert (ytl{1}, 'x1: 0.05 to 1');
%! assert (ytl{2}, 'x2: 0.05 to 20');
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
%! assert (numel (h), 4);
%! assert (xd1(1), 8.10687671732127, -1e-10);
%! assert (xd1(2), -25.4487243632125, -1e-10);
%! assert (xd1(3), 0.661302203942261, -1e-10);
%! assert (yd1, [1 2 3]);
%! assert (xd2(1), 0.565266595687836, -1e-10);
%! assert (xd2(2), 15.6484868389547, -1e-10);
%! assert (yd2, [1 1]);
%! assert (xd3(1), -33.3368582824351, -1e-10);
%! assert (xd3(2), -17.5605904439899, -1e-10);
%! assert (yd3, [2 2]);
%! assert (xd4(1), -1.25582490831999, -1e-10);
%! assert (xd4(2), 2.57842931620451, -1e-10);
%! assert (yd4, [3 3]);
%! assert (get (ax, 'YTick'), [1 2 3]);
%! assert (ytl{1}, 'x1: 0.05 to 1');
%! assert (ytl{2}, 'x2: 0.05 to 20');
%! assert (ytl{3}, 'x3: 1.22465e-16 to 1');
%! assert (mean (xd2), xd1(1), 1e-10);
%! assert (mean (xd3), xd1(2), 1e-10);
%! assert (mean (xd4), xd1(3), 1e-10);
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
%! assert (numel (h), 3);
%! assert (xd1(1), 2.50035744908398, -1e-10);
%! assert (xd1(2), -19.5912988214488, -1e-10);
%! assert (yd1, [1 2]);
%! assert (xd2(1), 1.40421088339552, -1e-10);
%! assert (xd2(2), 3.59650401477245, -1e-10);
%! assert (yd2, [1 1]);
%! assert (xd3(1), -20.6333076647782, -1e-10);
%! assert (xd3(2), -18.5492899781194, -1e-10);
%! assert (yd3, [2 2]);
%! assert (ytl{1}, 'x1: 0.05 to 1');
%! assert (ytl{2}, 'x2: 0.05 to 20');
%! assert (mean (xd2), xd1(1), 1e-10);
%! assert (mean (xd3), xd1(2), 1e-10);
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
%! assert (numel (h), 3);
%! assert (xd1(1), 2.51587141860715, -1e-10);
%! assert (xd1(2), -19.6411669663483, -1e-10);
%! assert (yd1, [1 2]);
%! assert (xd2(1), 1.08491557053384, -1e-10);
%! assert (xd2(2), 3.94682726668046, -1e-10);
%! assert (yd2, [1 1]);
%! assert (xd3(1), -20.8383905241664, -1e-10);
%! assert (xd3(2), -18.4439434085302, -1e-10);
%! assert (yd3, [2 2]);
%! assert (ytl{1}, 'x1: 0.05 to 1');
%! assert (ytl{2}, 'x2: 0.05 to 20');
%! assert (mean (xd2), xd1(1), 1e-10);
%! assert (mean (xd3), xd1(2), 1e-10);
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
%! assert (numel (h), 3);
%! assert (xd1(1), 2.81335053251731, -1e-10);
%! assert (xd1(2), -19.8951125513936, -1e-10);
%! assert (yd1, [1 2]);
%! assert (xd2(1), 2.36234515544211, -1e-10);
%! assert (xd2(2), 3.26435590959250, -1e-10);
%! assert (yd2, [1 1]);
%! assert (xd3(1), -20.4919734818287, -1e-10);
%! assert (xd3(2), -19.2982516209584, -1e-10);
%! assert (yd3, [2 2]);
%! assert (ytl{1}, 'x1: 0.05 to 1');
%! assert (ytl{2}, 'x2: 0.05 to 20');
%! assert (mean (xd2), xd1(1), 1e-10);
%! assert (mean (xd3), xd1(2), 1e-10);
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! ax = axes (fig);
%! h = plotEffects (ax, mdl);
%! assert (isequal (get (h(1), 'Parent'), ax));
%! assert (get (h(1), 'XData'), [2.38302891604232, -19.5277648300125], -1e-10);
%! close (fig);

%!test
%! fig = figure ('visible', 'off');
%! h = plotEffects (mdl);
%! assert (isequal (get (h(1), 'Parent'), gca ()));
%! assert (get (h(1), 'XData'), [2.38302891604232, -19.5277648300125], -1e-10);
%! close (fig);

%!test
%! ## numeric predictor adjusted data
%! fig = figure ('visible', 'off');
%! h = plotAdjustedResponse (mdl, 'x1');
%! assert (get (h(1), 'XData'), X(:,1)', 1e-10);
%! assert (get (h(1), 'YData'), ...
%!   [-6.70590752804287, -6.54551694016554, -6.5544435914624,  -6.59143572669715, ...
%!    -6.49138418414686, -6.21712299745016, -5.89359961368025, -5.69299878673044, ...
%!    -5.67643670865536, -5.73777106454765, -5.70118778736348, -5.48284370035446, ...
%!    -5.16795154710756, -4.93243578806991, -4.88118846293094, -4.95163193306634, ...
%!    -4.97125247389768, -4.81620859768203, -4.52519034621851, -4.26384784484646], 1e-8);
%! xf = get (h(2), 'XData');
%! yf = get (h(2), 'YData');
%! assert (numel (xf), 100);
%! assert (xf(1:5), ...
%!   [0.05, 0.0595959595959596, 0.0691919191919192, 0.0787878787878788, 0.0883838383838384], 1e-10);
%! assert (yf(1:5), ...
%!   [-6.78153223917696, -6.75746124002502, -6.73339024087308, -6.70931924172113, -6.68524824256919], 1e-8);
%! assert (xf(end-4:end), ...
%!   [0.961616161616162, 0.971212121212121, 0.980808080808081, 0.99040404040404, 1], 1e-10);
%! assert (yf(end-4:end), ...
%!   [-4.49478731974241, -4.47071632059047, -4.44664532143853, -4.42257432228658, -4.39850332313464], 1e-8);
%! close (fig);

%!test
%! ## title and axis labels follow the standard convention
%! fig = figure ('visible', 'off');
%! plotAdjustedResponse (mdl, 'x1');
%! assert (get (get (gca, 'Title'),  'String'), 'Adjusted Response Plot');
%! assert (get (get (gca, 'XLabel'), 'String'), 'x1');
%! assert (get (get (gca, 'YLabel'), 'String'), 'Adjusted y');
%! close (fig);

%!test
%! ## ax routing and a second predictor
%! fig = figure ('visible', 'off');
%! ax  = axes (fig);
%! h   = plotAdjustedResponse (ax, mdl, 'x2');
%! assert (isequal (get (h(1), 'Parent'), ax));
%! assert (get (h(1), 'XData'), X(:,2)', 1e-8);
%! assert (get (h(1), 'YData'), ...
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
%! assert (get (h(1), 'Marker'), 's');
%! assert (get (h(1), 'MarkerSize'), 10);
%! assert (get (h(1), 'Color'), [1 0 0]);
%! assert (get (h(2), 'Marker'), 'none');
%! close (fig);

%!test
%! yn = y;
%! yn(3) = NaN;
%! mn = fitlm (X, yn);
%! fig = figure ('visible', 'off');
%! h = plotAdjustedResponse (mn, 'x1');
%! xd = get (h(1), 'XData');
%! yd = get (h(1), 'YData');
%! assert (isnan (xd(3)));
%! assert (isnan (yd(3)));
%! assert (yd([1 2 4]), [-7.04679028889336, -6.88651148335619, -6.93287739924846], 1e-8);
%! close (fig);

%!test
%! w  = mod ((1:n)', 3) + 1;
%! mw = fitlm (X, y, 'Weights', w);
%! fig = figure ('visible', 'off');
%! h = plotAdjustedResponse (mw, 'x1');
%! assert (get (h(1), 'YData'), ...
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
%! assert (get (h(1), 'YData'), ...
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
%! assert (get (h(1), 'XData'), wt', 1e-10);
%! assert (get (h(1), 'YData'), ...
%!   [22.1247351073949, 19.1247351073949, 22.1247351073949, 20.1247351073949, ...
%!    21.1247351073949, 18.6102199722546, 20.6102199722546, 22.6102199722546, ...
%!    25.8864161365654, 27.8864161365654, 23.8864161365654, 21.8864161365654], 1e-8);
%! xf = get (h(2), 'XData');
%! yf = get (h(2), 'YData');
%! assert (xf(1:5), ...
%!   [2500, 2512.17171717172, 2524.34343434343, 2536.51515151515, 2548.68686868687], 1e-8);
%! assert (yf(1:5), ...
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
%! assert (get (h(1), 'XData'), [1 1 1 1 1 2 2 2 3 3 3 3]);
%! assert (get (h(1), 'YData'), ...
%!   [19.2479799272553, 17.3911214761304, 18.8366909043795, 16.8185458004291, ...
%!    17.9153196881646, 22.2641057484776, 24.463701891932,  23.9415324428265, ...
%!    28.7560523180671, 27.1754184718549, 24.3850920685482, 24.8044392619348], 1e-7);
%! assert (get (h(2), 'XData'), [1 2 3]);
%! assert (get (h(2), 'YData'), [18.0419315592718, 23.5564466944121, 26.2802505301012], 1e-8);
%! assert (get (gca, 'XTickLabel'), {'70'; '76'; '82'});
%! close (fig);

%!test
%! ## added variable plot for the whole model
%! fig = figure ('visible', 'off');
%! h = plotAdded (mdl);
%! assert (get (get (gca, 'Title'), 'String'), 'Added Variable Plot for Whole Model');
%! assert (get (get (gca, 'XLabel'), 'String'), 'Adjusted Whole Model');
%! assert (get (h(1), 'XData'), ...
%!   [0.0284033824838481, 0.0204548552857604, -0.0238455815942649, -0.104497928156226, ...
%!    -0.221502184400124, -0.374858350325959, -0.564566425933731, -0.79062641122344, ...
%!    -1.05303830619508, -1.35180211084867, -1.68691782518418, -2.05838544920164, ...
%!    -2.46620498290103, -2.91037642628236, -3.39089977934563, -3.90777504209083, ...
%!    -4.46100221451797, -5.05058129662704, -5.67651228841805, -6.338795189891], 1e-7);
%! assert (get (h(1), 'YData'), ...
%!   [0.268294196961578, 0.281859485365135, 0.0282240016119717, -0.351360499061586, ...
%!    -0.691784854932629, -0.955883099639786, -1.26860268025624, -1.80212835067532, ...
%!    -2.61757630295165, -3.60880422217788, -4.59999804131014, -5.50731458360009, ...
%!    -6.41596659263467, -7.50187852886103, -8.86994243196858, -10.457580663333, ...
%!    -12.0922794983759, -13.6501974493543, -15.1700245580674, -16.8174109498545], 1e-7);
%! assert (get (h(2), 'DisplayName'), 'Fit: y = 2.69267*x');
%! close (fig);

%!test
%! ## added variable plot for just the intercept term
%! fig = figure ('visible', 'off');
%! h = plotAdded (mdl, 1);
%! assert (get (get (gca, 'Title'), 'String'), 'Added Variable Plot for (Intercept)');
%! assert (get (h(1), 'YData'), ...
%!   [-5.41993222738086, -5.40485070721962, -5.55724508427077, -5.73586360329798, ...
%!    -5.77559710257835, -5.63927961575051, -5.45185858988763, -5.38551877888305, ...
%!    -5.50137637479139, -5.6932890627053,  -5.78544277558092, -5.6939943366699,  ...
%!    -5.50415648955918, -5.3918536946959,  -5.4619779917695,  -5.65195174215564, ...
%!    -5.78926122127593, -5.75006494138741, -5.57305294428921, -5.42387535532067], 1e-7);
%! assert (get (h(2), 'DisplayName'), 'Fit: y = 0.116189*x');
%! close (fig);

%!test
%! ## added variable plot for one predictor picked by index
%! fig = figure ('visible', 'off');
%! h = plotAdded (mdl, 2);
%! assert (get (get (gca, 'Title'), 'String'), 'Added Variable Plot for x1');
%! assert (get (h(1), 'YData'), ...
%!   [-5.90289714234183, -5.75941203626873, -5.79651449057265, -5.87295275001727, ...
%!    -5.82361765287968, -5.6113432327985,  -5.36107693684693, -5.24500351891828, ...
%!    -5.32423917106718, -5.49264157838629, -5.57439667383173, -5.48566128065516, ...
%!    -5.31164814244354, -5.22828171964398, -5.34045405194592, -5.58558750072505, ...
%!    -5.79116834140295, -5.83335508623668, -5.75083777702536, -5.70926653910833], 1e-7);
%! assert (get (h(2), 'DisplayName'), 'Fit: y = 2.50845*x');
%! xf = get (h(2), 'XData');
%! yf = get (h(2), 'YData');
%! assert (xf(1:3), [0.370121951219512, 0.372449462532903, 0.374776973846294], 1e-9);
%! assert (yf(1:3), [-5.97852185347592, -5.97268340425253, -5.96684495502913], 1e-7);
%! close (fig);

%!test
%! ## added variable plot for the other predictor, a negative slope this time
%! fig = figure ('visible', 'off');
%! h = plotAdded (mdl, 3);
%! assert (get (get (gca, 'Title'), 'String'), 'Added Variable Plot for x2');
%! assert (get (h(1), 'YData'), ...
%!   [-8.3040737600235,  -7.38815394983204, -6.7394349117973,  -6.21666489068295, ...
%!    -5.65473472476609, -5.01647844768535, -4.4268435065139,  -4.05801465514508, ...
%!    -3.9711080856335,  -4.05998148307182, -4.14882078041619, -4.15378280091823, ...
%!    -4.16008028816491, -4.34363770260337, -4.80934708392301, -5.49463079349955, ...
%!    -6.22697510675454, -6.88253853594507, -7.50001112287024, -8.2450429928694], 1e-7);
%! assert (get (h(2), 'DisplayName'), 'Fit: y = -0.978835*x');
%! close (fig);

%!test
%! ## ax argument sends the plot to that axes
%! fig = figure ('visible', 'off');
%! ax  = axes (fig);
%! h   = plotAdded (ax, mdl, 2);
%! assert (isequal (get (h(1), 'Parent'), ax));
%! close (fig);

%!test
%! ## name-value styling only changes the data points, not the fit line
%! fig = figure ('visible', 'off');
%! h = plotAdded (mdl, 2, 'Marker', 's', 'MarkerSize', 10, 'Color', 'r');
%! assert (get (h(1), 'Marker'), 's');
%! assert (get (h(1), 'MarkerSize'), 10);
%! assert (get (h(1), 'Color'), [1 0 0]);
%! assert (get (h(2), 'Marker'), 'none');
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
%! assert (isnan (xd(3)));
%! assert (isnan (yd(3)));
%! close (fig);

%!test
%! ## weighted fit still gives a proper added variable plot
%! w  = mod ((1:n)', 3) + 1;
%! mw = fitlm (X, y, 'Weights', w);
%! fig = figure ('visible', 'off');
%! h = plotAdded (mw, 2);
%! assert (get (h(1), 'YData'), ...
%!   [-6.04899282212585, -5.90562605001686, -5.9429257275943,  -6.01964009962184, ...
%!    -5.97066000437658, -5.75881947549715, -5.50906596005672, -5.39358421194863, ...
%!    -5.4734904232275,  -5.64264227898597, -5.7252257121802,  -5.63739754606181, ...
%!    -5.46437052421778, -5.38206910709522, -5.49538533438357, -5.74174156745852, ...
%!    -5.94862408174164, -5.99219138949,    -5.91113353250271, -5.87110063611913], 1e-7);
%! assert (get (h(2), 'DisplayName'), 'Fit: y = 2.38836*x');
%! close (fig);

%!test
%! ## robust fit still gives a proper added variable plot
%! mr = fitlm (X, y, 'RobustOpts', 'on');
%! fig = figure ('visible', 'off');
%! h = plotAdded (mr, 2);
%! assert (get (h(1), 'YData'), ...
%!   [-5.89409568732029, -5.75060102429134, -5.78768755033552, -5.86410351021651, ...
%!    -5.81473974221138, -5.60243027995878, -5.35212257053188, -5.23600136782402, ...
%!    -5.31518286388981, -5.4835247438219,  -5.56521294057644, -5.47640427740507, ...
%!    -5.30231149789475, -5.2188590624926,  -5.33093901088806, -5.5759737044568,  ...
%!    -5.78144941862042, -5.82352466563597, -5.74088948730258, -5.69919400895959], 1e-7);
%! assert (get (h(2), 'DisplayName'), 'Fit: y = 2.48815*x');
%! close (fig);

%!test
%! ## added variable plot for the whole model
%! load carsmall
%! Year = categorical (Model_Year);
%! tbl  = table (MPG, Weight, Year);
%! mdl1  = fitlm (tbl, 'MPG ~ Year + Weight^2');
%! fig  = figure ('visible', 'off');
%! h    = plotAdded (mdl1);
%! assert (get (get (gca, 'Title'),  'String'), 'Added Variable Plot for Whole Model');
%! assert (get (get (gca, 'XLabel'), 'String'), 'Adjusted Whole Model');
%! xd = get (h(1), 'XData');
%! yd = get (h(1), 'YData');
%! assert (xd(1:5), ...
%!   [-4.54006581304461, -4.65629283432076, -4.49502736854252, ...
%!    -4.49300111649177, -4.50376945388336], 1e-7);
%! assert (yd(1:5), [18, 15, 18, 16, 17], 1e-8);
%! assert (get (h(2), 'DisplayName'), 'Fit: y = 8.44866*x');
%! xf = get (h(2), 'XData');
%! yf = get (h(2), 'YData');
%! assert (xf(1:3), [-5.06005142297709, -5.03050027197254, -5.000949120968], 1e-7);
%! assert (yf(1:3), [11.4556384009199, 11.7053059986795, 11.9549735964392], 1e-7);
%! close (fig);

%!test
%! ## added variable plot for the weight terms picked as a pair
%! load carsmall
%! Year = categorical (Model_Year);
%! tbl  = table (MPG, Weight, Year);
%! mdl1  = fitlm (tbl, 'MPG ~ Year + Weight^2');
%! fig  = figure ('visible', 'off');
%! h    = plotAdded (mdl1, [2 5]);
%! assert (get (get (gca, 'Title'),  'String'), 'Added Variable Plot for Specified Terms');
%! assert (get (get (gca, 'XLabel'), 'String'), 'Adjusted Specified Terms');
%! xd = get (h(1), 'XData');
%! yd = get (h(1), 'YData');
%! assert (xd(1:5), ...
%!   [-2181.48546848357, -2241.34798193195, -2158.28850051435, ...
%!    -2157.24488312395, -2162.79109548355], 1e-6);
%! assert (yd(1:5), ...
%!   [24.0284299339692, 21.0284299339692, 24.0284299339692, ...
%!    22.0284299339692, 23.0284299339692], 1e-7);
%! assert (get (h(2), 'DisplayName'), 'Fit: y = 0.0164036*x');
%! close (fig);

%!test
%! load hald
%! Xh = ingredients;
%! yh = heat;
%! mdlh = fitlm (Xh, yh);
%! assert (mdlh.Coefficients.Estimate, ...
%!   [62.405369299918; 1.55110264750845; 0.510167579684912; ...
%!    0.101909403579662; -0.144061029071018], 1e-9);
%! assert (mdlh.Coefficients.SE, ...
%!   [70.0709592085362; 0.744769867130993; 0.72378800183518; ...
%!    0.754709045051309; 0.70905206344651], 1e-9);
%! assert (mdlh.Coefficients.tStat, ...
%!   [0.890602469336764; 2.0826603169159; 0.704857746178952; ...
%!    0.135031379639465; -0.203174120065001], 1e-8);
%! assert (mdlh.Coefficients.pValue, ...
%!   [0.399133563385561; 0.0708216874297252; 0.500901103474289; ...
%!    0.895922690510107; 0.844071473291884], 1e-8);
%! assert (mdlh.CoefficientNames, {'(Intercept)', 'x1', 'x2', 'x3', 'x4'});
%! assert (mdlh.NumCoefficients,          5);
%! assert (mdlh.NumEstimatedCoefficients, 5);
%! assert (mdlh.DFE,                      8);
%! assert (mdlh.SSE, 47.863639350499,   1e-8);
%! assert (mdlh.SSR, 2667.89943757258,  1e-6);
%! assert (mdlh.SST, 2715.76307692308,  1e-6);
%! assert (mdlh.MSE,  5.98295491881254, 1e-9);
%! assert (mdlh.RMSE, 2.44600795559061, 1e-9);
%! assert (mdlh.Rsquared.Ordinary, 0.98237562040768, 1e-9);
%! assert (mdlh.Rsquared.Adjusted, 0.97356343061152, 1e-9);
%! assert (mdlh.LogLikelihood, -26.918344895826, 1e-8);
%! assert (mdlh.ModelCriterion.AIC,  63.8366897916521, 1e-7);
%! assert (mdlh.ModelCriterion.AICc, 72.4081183630806, 1e-7);
%! assert (mdlh.ModelCriterion.BIC,  66.6614365789598, 1e-7);
%! assert (mdlh.ModelCriterion.CAIC, 71.6614365789598, 1e-7);
%! assert (mdlh.Fitted(1:5), ...
%!   [78.4952395815018; 72.7887993002909; 105.970937532083; ...
%!    89.3271002550427; 95.649244438227], 1e-8);
%! assert (mdlh.Residuals.Raw(1:5), ...
%!   [0.00476041849822195; 1.51120069970906; -1.67093753208295; ...
%!    -1.72710025504266; 0.250755561773033], 1e-8);
%! assert (mdlh.Residuals.Pearson(1:5), ...
%!   [0.00194619910672879; 0.617823297040002; -0.683128412670876; ...
%!    -0.706089385807266; 0.102516249466771], 1e-8);
%! assert (mdlh.Residuals.Studentized(1:5), ...
%!   [0.00271470565323249; 0.734526653667679; -1.05809320265782; ...
%!    -0.824036396702643; 0.119767490249399], 1e-8);
%! assert (mdlh.Residuals.Standardized(1:5), ...
%!   [0.00290214088954622; 0.756624558354514; -1.05027405557414; ...
%!    -0.841081414787206; 0.127905848829164], 1e-8);
%! assert (mdlh.Diagnostics.Leverage(1:5), ...
%!   [0.550284813713987; 0.333242829857405; 0.576942476415795; ...
%!    0.29523667959374; 0.357601364034465], 1e-8);
%! assert (mdlh.Diagnostics.CooksDistance(1:5), ...
%!   [2.06118491039641e-06; 0.0572247602223712; 0.300862709270433; ...
%!    0.0592697490074745; 0.00182140011900327], 1e-8);
%! assert (mdlh.Diagnostics.Dffits(1:5), ...
%!   [0.00300294746488125; 0.519283016757055; -1.23563576459509; ...
%!    -0.533347058289184; 0.0893585735072963], 1e-7);
%! assert (mdlh.Diagnostics.S2_i(1:5), ...
%!   [6.83765556564705; 6.34835899957971; 5.89485540180634; ...
%!    6.23302709557492; 6.82367982420565], 1e-7);
%! assert (mdlh.Diagnostics.CovRatio(1:5), ...
%!   [4.33530738335252; 2.01725612858557; 2.19476339013102; ...
%!    1.74129811023362; 3.00406926806094], 1e-6);
%! assert (coefCI (mdlh), ...
%!   [-99.178552392689 223.989290992525; -0.166339745871082 3.26854504088797; ...
%!    -1.15889054555817 2.179225704928; -1.63845277518465 1.84227158234397; ...
%!    -1.77913801945372 1.49101596131168], 1e-7);
%! assert (coefCI (mdlh, 0.1), ...
%!   [-67.8949453842232 192.705683984059; 0.166167302672858 2.93603799234403; ...
%!    -0.835750978716108 1.85608613808593; -1.30150832005232 1.50532712721164; ...
%!    -1.46257740216021 1.17445534401817], 1e-7);
%! [p, F, r] = coefTest (mdlh);
%! assert (p, 4.75618174559791e-07, 1e-12);
%! assert (F, 111.479171821258, 1e-6);
%! assert (r, 4);
%! [dw, pdw] = dwtest (mdlh);
%! assert (dw,  0.842123108585363, 1e-9);
%! assert (pdw, 2.05259693286049,  1e-8);

%!test
%! load hald
%! Xq = ingredients;
%! yq = heat;
%! mdlq = fitlm (Xq, yq, 'purequadratic');
%! assert (mdlq.Coefficients.Estimate, ...
%!   [-210.864812527187; 4.01775196981369; 5.27927179495849; 4.98703005469684; ...
%!    1.18967556414545; -0.00475259718542981; -0.0278555063988026; ...
%!    -0.0885225308739459; 0.0125632231921875], 1e-8);
%! assert (mdlq.Coefficients.SE, ...
%!   [62.2492955088454; 0.709968572684776; 0.928573133899706; 1.15325562045609; ...
%!    0.559753134550446; 0.0119258668009245; 0.00570193167405062; ...
%!    0.0224063172978396; 0.0040139580330777], 1e-8);
%! assert (mdlq.Coefficients.tStat, ...
%!   [-3.38742488253901; 5.6590560827508; 5.68535918413584; 4.3243058748108; ...
%!    2.12535757410437; -0.398511677579811; -4.88527537528597; ...
%!    -3.95078449069724; 3.12988404180068], 1e-7);
%! assert (mdlq.Coefficients.pValue, ...
%!   [0.0275955163014823; 0.00480588842249351; 0.00472567672152044; ...
%!    0.0124052274432915; 0.100732674864514; 0.710609610173011; ...
%!    0.00812965117162355; 0.0168069450605605; 0.0351891921636521], 1e-7);
%! assert (mdlq.CoefficientNames, ...
%!   {'(Intercept)', 'x1', 'x2', 'x3', 'x4', 'x1^2', 'x2^2', 'x3^2', 'x4^2'});
%! assert (mdlq.NumCoefficients, 9);
%! assert (mdlq.DFE, 4);
%! assert (mdlq.Rsquared.Ordinary, 0.998060528051984, 1e-9);
%! assert (mdlq.Rsquared.Adjusted, 0.994181584155951, 1e-9);
%! assert (mdlq.SSE, 5.26714630515049, 1e-7);
%! assert (mdlq.SSR, 2710.49593061793, 1e-6);
%! assert (mdlq.SST, 2715.76307692308, 1e-6);
%! Xnewq = mean (Xq, 1);
%! [ypredq, yciq] = predict (mdlq, Xnewq);
%! assert (ypredq, 101.90428629036, 1e-7);
%! assert (yciq, [97.996264336243, 105.812308244477], 1e-6);
%! [ypredq2, yciq2] = predict (mdlq, Xnewq, 'Alpha', 0.01);
%! assert (ypredq2, 101.90428629036, 1e-7);
%! assert (yciq2, [95.4237317847484, 108.384840795971], 1e-6);
%! yfeq = feval (mdlq, Xnewq(1), Xnewq(2), Xnewq(3), Xnewq(4));
%! assert (yfeq, 101.90428629036, 1e-7);
%! ysimq = random (mdlq, Xnewq);
%! assert (isscalar (ysimq));
%! assert (isnumeric (ysimq));
%! mdlq2 = removeTerms (mdlq, 'x1^2');
%! assert (mdlq2.Coefficients.Estimate, ...
%!   [180.815670525319; 0.228012666917658; -2.0654020024496; -1.90840414992816; ...
%!    -0.0108001957000414; 0.0257318558607856; 0.00699693273966028], 1e-8);
%! assert (mdlq2.CoefficientNames, ...
%!   {'(Intercept)', 'x2', 'x3', 'x4', 'x2^2', 'x3^2', 'x4^2'});
%! assert (mdlq2.NumCoefficients, 7);
%! mdlq3 = addTerms (mdlq2, 'x1^2');
%! assert (mdlq3.Coefficients.Estimate, mdlq.Coefficients.Estimate, 1e-8);
%! assert (mdlq3.CoefficientNames, mdlq.CoefficientNames);
%! assert (mdlq3.SSE, 5.26714630515049, 1e-7);

%!test
%! load hald
%! Xr = ingredients;
%! yr = heat;
%! mdlr = fitlm (Xr, yr, 'RobustOpts', 'bisquare');
%! assert (mdlr.Coefficients.Estimate, ...
%!   [60.0897358816096; 1.57529551556915; 0.532199192097796; ...
%!    0.133455378556458; -0.120521170556001], 1e-8);
%! assert (mdlr.Coefficients.SE, ...
%!   [75.8175597390933; 0.805849306629754; 0.783146694256936; ...
%!    0.816603608044244; 0.767202244491812], 1e-8);
%! assert (mdlr.Coefficients.tStat, ...
%!   [0.792556976093573; 1.95482642053437; 0.679565138945976; ...
%!    0.163427368238162; -0.157091785668371], 1e-7);
%! assert (mdlr.Coefficients.pValue, ...
%!   [0.450897370203866; 0.0863457969332376; 0.515957116726031; ...
%!    0.874235088124976; 0.879064839096153], 1e-7);
%! assert (mdlr.CoefficientNames, {'(Intercept)', 'x1', 'x2', 'x3', 'x4'});
%! assert (is_function_handle (mdlr.Robust.RobustWgtFun));
%! assert (mdlr.Robust.Tune, 4.685, 1e-10);
%! assert (size (mdlr.Robust.Weights), [13, 1]);
%! assert (isnumeric (mdlr.Robust.Weights));
%! assert (mdlr.SSE, 56.0362670671825, 1e-6);
%! assert (mdlr.MSE, 7.00453338339782, 1e-8);
%! assert (mdlr.RMSE, 2.64660790133291, 1e-8);
%! assert (mdlr.Rsquared.Ordinary, 0.97929734395902, 1e-9);
%! assert (mdlr.Rsquared.Adjusted, 0.96894601593853, 1e-9);
%! assert (mdlr.DFE, 8);
%! H = [0 1 -1 0 0];
%! [p1, F1, r1] = coefTest (mdlr, H);
%! assert (p1, 0.00308748318894346, 1e-11);
%! assert (F1, 17.4568343157849, 1e-7);
%! assert (r1, 1);
%! [pd1, dw1] = dwtest (mdlr, 'exact', 'both');
%! assert (pd1, 0.844119247360191, 1e-9);
%! assert (dw1, 2.05387711905232, 1e-8);
%! [pd2, dw2] = dwtest (mdlr, 'approximate', 'right');
%! assert (pd2, 0.425180546504485, 1e-9);
%! assert (dw2, 2.05387711905232, 1e-8);
%! Xnewr = mean (Xr, 1);
%! [ypredr, ycir] = predict (mdlr, Xnewr, 'Simultaneous', true);
%! assert (ypredr, 95.4263340097424, 1e-7);
%! assert (ycir, [92.2744598719279, 98.5782081475569], 1e-6);
%! [ypredr2, ycir2] = predict (mdlr, Xnewr, 'Simultaneous', true, 'Alpha', 0.1);
%! assert (ypredr2, 95.4263340097424, 1e-7);
%! assert (ycir2, [92.7161333049321, 98.1365347145527], 1e-6);

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
%!error <indexing is not supported> mdl (1)
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
