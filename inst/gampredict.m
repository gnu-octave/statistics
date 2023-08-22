## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
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

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{yFit} =} gampredict (@var{X}, @var{Y}, @var{Xfit})
## @deftypefnx {statistics} {@var{yFit} =} gampredict (@var{X}, @var{Y}, @var{Xfit}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{yFit}, @var{ySD}, @var{yInt}] =} gampredict (@var{X}, @var{Y}, @var{Xfit}, @var{name}, @var{value})
##
## Predict new data points using trained GAM model, by GAM regression
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the GAM model.
## @item
## @code{Y} is @math{Nx1} numeric matrix containing the Response data for
## corresponding predictor data in @var{X}.
## @var{Y} must have same numbers of Rows as @var{X}.
## @item
## @code{Xfit} must be a @math{MxP} numeric matrix of query/new points to
## predict the response value.
## @var{Xfit} must have same numbers of columns as @var{X}.
##
## @emph{Additional parameters can be passed as name value pairs :}
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @tab @qcode{"formula"} @tab  model specification specified as a
## string of the form @qcode{'Y ~ terms'} where 'Y' represents the reponse
## variable and 'terms' the predictor variables. formula is used to specify
## the subset of variables for training model.
## for example @qcode{"Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3"} specifies the
## linear terms x1, x2, x3, and x4 for predictor variables, x1:x2 and x2:x3
## specifies the interaction term for x1, x2 and x2, x3 respectively.
##
## @item @tab @qcode{"responsename"} @tab Response Variable Name specified as
## a string. default value is 'Y'.
##
## @item @tab @qcode{"predictors"} @tab Predictor Variable names, specified as
## cell of string(s). the length or columns of @qcode{"predictors"} must be
## same as @qcode{"X"}. If not supplied the program will generate default
## variable names (x1, x2, ... xn) for each column in @qcode{"X"}.
##
## @item @tab @qcode{"fitstd"} @tab Logical Value 0(false) or 1(true) to
## specify flag to fit model for the standard deviation of the response
## variable.
##
## @item @tab @qcode{"interactions"} @tab
##
## @item @tab @qcode{"maxpval"} @tab Maximum p-value for detecting interaction
## terms, must be a numeric scalar between 0 to 1. Interaction terms with
## p-vale less than maxpval will be used. default value is set to 0.05.
##
## @item @tab @qcode{"catpredictors"} @tab List of categorical predictors in
## predictor data @qcode{"X"}, specified as the index of column in @qcode{"X"}.
##
## @item @tab @qcode{"Weights"} @tab Observational Weights specified as a
## numeric matrix with each row correspoding to the observations in @qcode{"X"}.
## @qcode{"Weights"} must have same number of rows as @qcode{"X"}. Default is
## ones (size (X,1),1).
##
## @item @tab @qcode{"Alpha"} @tab Significance level of the prediction
## intervals @qcode{"yInt"}. Specified as scalar in range [0,1]. This argument
## is only valid when @qcode{"fitstd"} is set true. default value is 0.05.
## for example 'Alpha',0.05 return 95% prediction intervals.
##
## @item @tab @qcode{"dof"} @tab Degree of freedom to fit a third order spline.
## for fitting a spline @qcode{"dof = knots + order"}, for fitting a GAM a
## polynomial spline of degree '3' is used hence the number of knots can be
## controlled by degree of freedom, degree of freedom can be used to adjust the
## fit of the each variable. the length of @qcode{"dof"} must be same as the
## columns of @qcode{"X"}. default value is 8 for each predictor variable.
##
## @item @tab @qcode{"includeinteractions"} @tab flag to indicate weather to
## include interactions to predict new values from.
##
## @end multitable
##
## @code{@var{yFit} = gampredict (@var{X}, @var{Y}, @var{Xfit})} returns a
## Numeric Matrix of predicted Response values for predictor data in @var{Xfit},
## using a generalised additive model trained using predictor values in @var{X}
## and response values in @var{Y}.
##
##
## @code{@var{yFit} = gampredict (@var{X}, @var{Y}, @var{Xfit}, @var{name}, @var{value})}
## returns a matrix of predicted Response values for predictor data
## in @var{Xfit} using GAM regression model with specified options as Name-value
## pairs.
##
## @code{[@var{yFit}, @var{ySD}, @var{yInt}] = gampredict (@var{X}, @var{Y}, @var{Xfit})}
## returns @var{yFit} containing the predicted response values for predictor
## data in @var{Xfit}, It also returns :
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"ySD"} @tab Numeric matrix of Standard deviations of
## predicted reponse variable for each observation.
##
## @item @tab @qcode{"yInt"} @tab Prediction intervals of response variable as
## @math{nx2} numeric matrix with first column containing the lower limits and
## second column the upper limit of prediction intervals.
## Each row of the prediction intervals contains @math{100 x (1-Alpha) %}
## of prediction interval.
##
##@end multitable
##
## @end itemize
## for demo use demo gampredict
##
## @seealso{regress}
## @end deftypefn

function [yFit, varargout] = gampredict (X, Y, Xfit, varargin)

  ## Check minimum number of input arguments
  if (nargin < 3)
	  error ("gampredict: too few input arguments.");
  endif

  ## Get training sample size and number of variables in training data
  nsample = rows (X);
  ndims_X = columns (X);

  ## Check size consistency in data
  if (nsample != rows (Y))
    error ("gampredict: number of rows in X and Y must be equal.");
  endif

  if (ndims_X != columns (Xfit))
    error ("gampredict: number of columns in Xfit must be equal to X.");
  endif

  ## Process optional parameters

  ## Default values for optional parameters
  Alpha = 0.05;                   # significance level for intervals
  Formula = [];                   # formula for GAM model
  Interactions = [];               # number or list of interaction terms
  MaxPValue = 1;                  # max p-value for including interaction terms
  Knots = ones (1, ndims_X) * 5;  # Knots
  DoF = ones (1, ndims_X) * 8;    # Degrees of freedom for fitting spline
  Order = ones (1, ndims_X) * 3;  # Order of spline
  Tol = 1e-2;                     # tolerence for converging splines
  PredictorNames = [];
  ResponseName   = [];

  ## Number of arguments for Knots, DoF, Order (maximum 2 allowed)
  KDO = 0;

  while (numel (varargin) > 0)
    switch (tolower (varargin {1}))
      case "alpha"
        Alpha = varargin{2};
        if (! (isnumeric (Alpha) && isscalar (Alpha) && Alpha > 0 && Alpha < 1))
          error ("gampredict: alpha must be numeric and in between 0 and 1.");
        endif

      case "formula"
        Formula = varargin{2};
        if (! ischar (Formula))
          error ("gampredict: formula must be a string.");
        endif

      case "interactions"
        tmp = varargin{2};
        if (isnumeric (tmp) && isscalar (tmp) && tmp == fix (tmp) && tmp >= 0)
          Interactions = tmp;
        elseif (ismatrix (tmp) && islogical (tmp))
          Interactions = tmp;
        elseif (ischar (tmp) && strcmpi (tmp, "all"))
          Interactions = tmp;
        else
          error ("gampredict: wrong value for interactions parameter.");
        endif

      case "maxpval"
        MaxPValue = varargin{2};
        if (! isnumeric (MaxPValue) ||
            ! isscalar (MaxPValue) || MaxPValue < 0 || MaxPValue > 1)
          error ("gampredict: MaxPValue must range from 0 to 1.");
        endif

      case "knots"
        if (KDO < 2)
          Knots = varargin{2};
          if (! isnumeric (Knots) ||
              ! (isscalar (Knots) || isequal (size (Knots), [1, ndims_X])))
            error ("gampredict: invalid value for Knots.");
          endif
          DoF = Knots + Order;
          Order = DoF - Knots;
          KDO += 1;
        else
          error ("gampredict: DoF and Order have been set already.");
        endif

      case "dof"
        if (KDO < 2)
          DoF = varargin{2};
          if (! isnumeric (DoF) ||
              ! (isscalar (DoF) || isequal (size (DoF), [1, ndims_X])))
            error ("gampredict: invalid value for DoF.");
          endif
          Knots = DoF - Order;
          Order = DoF - Knots;
          KDO += 1;
        else
          error ("gampredict: Knots and Order have been set already.");
        endif

      case "order"
        if (KDO < 2)
          Order = varargin{2};
          if (! isnumeric (Order) ||
              ! (isscalar (Order) || isequal (size (Order), [1, ndims_X])))
            error ("gampredict: invalid value for Order.");
          endif
          DoF = Knots + Order;
          Knots = DoF - Order;
          KDO += 1;
        else
          error ("gampredict: DoF and Knots have been set already.");
        endif

      case "tol"
        Tol = varargin{2};
        if (! isnumeric (Tol) || ! isscalar (Tol) || !(Tol > 0))
          error ("gampredict: Tolerence (tol) must be a scalar > 0.");
        endif
      case "responsename"
        ResponseName = varargin{2};
        if (! ischar (ResponseName))
          error ("gampredict: ResponseName must be a string or char.");
        endif
      case "predictors"
        PredictorNames = varargin{2};
        if (! isempty (PredictorNames))
          if (! iscellstr (PredictorNames))
            error ("gampredict: PredictorNames must be supplied as cellstr.");
          elseif (columns (PredictorNames) != columns (X))
            error ("gampredict: PredictorNames must have same number of columns as X.");
          endif
        endif

      otherwise
        error ("gampredict: invalid NAME in optional pairs of arguments.")
    endswitch
    varargin (1:2) = [];
  endwhile


  ## Process predictor matrix
  notnans  = ! logical (sum (isnan ([Y, X]), 2));
  notnansf = ! logical (sum (isnan (Xfit), 2));
  Y        = Y (notnans);
  X        = X (notnans, :);
  Xfit     = Xfit (notnansf, :);


  ## Adding interaction terms to the predictor matrix for training
  if (isempty (PredictorNames))
    ## empty predictornames generate default
    PredictorNames = gendefault (columns (X));
  endif
  IntMat = [];
  if (! isempty (Formula))
    IntMat = parseInteractions (formula, PredictorNames);
    if (isempty (intMat))
      ## user has not Specified any Ineractions in the Formula
      ## check if specified in interactions
      if (! isempty (Interactions))
        IntMat = parseInteractions (Interactions, PredictorNames);
      endif
    endif
  else
    ## formula is empty train model for given terms and check Interactions
    if (! isempty (Interactions))
      IntMat = parseInteractions (Interactions, PredictorNames);
    endif
  endif

  if (! isempty (IntMat))
    for i = 1:rows (IntMatat)
      Xint = X (:,IntMat (i, 1)) .* X (:,IntMat (i, 1));
      X    = [X, Xint]; ## adding interaction term column
    endfor
  endif


  ## Fit the GAM model

  ## Compute intercept and residuals
  intercept   = mean (Y);
  res         = Y - intercept;
  converged   = false;
  num_itn     = 0;
  RSS         = zeros (1, columns (X));

  ## Main loop
  while (! converged)
    num_itn += 1;
    ## Single cycle of backfit
    for j = 1:columns (X)

      ## Calculate residuals to fit spline
      if (num_itn > 1)
        res = res + ppval (ppfk (j), X(:, j));
      endif

      ## Fit an spline to the data
      ## You used `order(j)-1` previously. Why?
      gk = splinefit (X(:,j), res, Knots(j), 'order', Order(j));

      ## This might be wrong! We need to check this out
      RSSk (j) = abs (sum (abs (Y - ppval (gk, X(:,j)) ...
                                  - intercept )) .^ 2) / nsample;
      ppfk (j) = gk;
      res = res - ppval (ppfk (j), X (:,j));
    endfor

    ## check if RSS is less than the tolerence
    if (all (abs (RSS - RSSk) <= Tol))
      converged = true;
    endif

    ## update RSS
    RSS = RSSk;
  #until (converged)
  endwhile


  ## Predict values from testing data
  yFit = predict_val (ppfk, Xfit, intercept);

  ## Predict Standard Deviation and Intervals of estimated data (if requested)
  if (nargout > 0)
    ## Predict response from training predictor data with the trained model
    yrs = predict_val (ppfk, X , intercept);
    ## Get the residuals between predicted and actual response data

    ## You need to check this code here !!
    rs     = Y - yrs;
    var_rs = var (rs);    ##
    var_pr = var (yFit);  ## var is calculated here instead take sqrt(SD)

    t_mul  = tinv (1 - Alpha / 2, DoF);


    ydev   = (yFit - mean (yFit)) .^ 2;
    ySD    = sqrt (ydev / (rows (yFit) - 1));

    varargout{1} = ySD;
    if (nargout > 1)
      #moe    = t_mul (1) * sqrt (var_pr + var_rs)
      moe    = t_mul (1) * ySD;
      lower  = (yFit - moe);
      upper  = (yFit + moe);
      yInt   = [lower, upper];
      varargout{2} = yInt;
    endif
  endif

endfunction

## Helper function for making prediction of new data based on GAM model
function ypred = predict_val (ppfk, X, intercept)
  [nsample, ndims_X] = size (X);
  ypred = ones (nsample, 1) * intercept;
  ## Add the remaining terms
  for j = 1:ndims_X
    ypred = ypred + ppval (ppfk(j), X (:,j));
  endfor
endfunction

## Function parseInteractions to detect user given interaction
##
## This function has some logic errors and needs updating!!
## It's not being used at the moment, so I leave it here as is
##
function intMat = parseInteractions (Formula, PredictorNames)
  intMat = [];
  if (islogical (Formula))
    if (numel (PredictorNames) != columns (Formula))
      error ("gampredict: ");
    endif
    for i = 1:rows (Formula)
      ind = find (Formula(i, :) == 1);
      intMat = [intMat; ind];
    endfor
  elseif (ischar (Formula))
    if (strcmpi (Formula, "all"))
      #calculate all p*(p-1)/2 interaction terms

    else
      #calculate interaction matrix by formula
      formulaParts = strsplit(Formula, '~');
      responseVar = strtrim(formulaParts{1});
      predictorString = strtrim(formulaParts{2});
      intterms = strsplit(predictorString, '+');
      in = [];
      for i = 1:numel (intterms)
        if (cell2mat (strfind (intterms (i), ':')))
          inters = strtrim (strsplit (cell2mat (intterms (i)), ':'));
          in = [in; (index (PredictorNames, inters(1)) + ...
                                          index (PredictorNames, inters(2)))];
        endif
      endfor
      ## pass recursively to get intMat
      intMat = parseInteractions (logical(in), PredictorNames);
    endif
  else
    error ("gampredict: Invalid value in Interactions.");
  endif
endfunction


##------ gendefault to generate default predictornames-----##
function defs = gendefault (p)
  for i = 1:p
    defs {i} = strcat ('x', num2str (i));
  endfor
endfunction


%!demo
%! ## Declare two different functions
%! f1 = @(x) cos (3 * x);
%! f2 = @(x) x .^ 3;
%!
%! ## Generate 80 samples for f1 and f2
%! x = [-4*pi:0.1*pi:4*pi-0.1*pi]';
%! X1 = f1 (x);
%! X2 = f2 (x);
%!
%! ## Create a synthetic response by adding noise
%! rand ("seed", 3);
%! Ytrue = X1 + X2;
%! Y = Ytrue + Ytrue .* 0.2 .* rand (80,1);
%!
%! ## Assemble predictor data
%! X = [X1, X2];
%!
%! ## Train the GAM and test on the same data
%! [ypred, ySDsd, yInt] = gampredict (X, Y, X, "order", [5,5]);
%!
%! ## Plot the results
%! figure
%! [sortedY, indY] = sort (Ytrue);
%! plot (sortedY, 'r-');
%! xlim ([0, 80]);
%! hold on
%! plot (ypred(indY), "g+")
%! plot (yInt(indY,1),'k:')
%! plot (yInt(indY,2),'k:')
%! xlabel ("Predictor samples");
%! ylabel ("Response");
%! title ("actual vs predicted values for function f1(x) = cos (3x) ");
%! legend ({"Theoretical Response", "Predicted Response", "Prediction Intervals"});
%!
%! ## Use 30% Holdout partitioning for training and testing data
%! C = cvpartition (80, "HoldOut", 0.3);
%! [ypred, ySDsd, yInt] = gampredict (X(training(C),:), Y(training(C)), X(test(C),:));
%!
%! ## Plot the results
%! figure
%! [sortedY, indY] = sort (Ytrue(test(C)));
%! plot (sortedY, 'r-');
%! xlim ([0, sum(test(C))]);
%! hold on
%! plot (ypred(indY), "g+")
%! plot (yInt(indY,1),'k:')
%! plot (yInt(indY,2),'k:')
%! xlabel ("Predictor samples");
%! ylabel ("Response");
%! title ("actual vs predicted values for function f1(x) = cos (3x) ");
%! legend ({"Theoretical Response", "Predicted Response", "Prediction Intervals"});


%!demo
%! ## Generate 101 samples based on a sinusoidal function and add some noise
%! X = [0:0.1*pi:2*pi-0.1*pi]';
%! Y = sin (X);
%! rand ("seed", 3);
%! Y += Y .* 0.5 .* (rand (20, 1) - 0.5);
%!
%! ## Train the GAM and test on the same data
%! [ypred, ySDsd, yInt] = gampredict (X, Y, X);
%! [Ysorted, I] = sort (Y);
%!
%! ## Plot the results
%! figure
%! plot (X, ypred, '*', X, Y, '+', X, sin (X), '-')
%! ylim ([-1.2, 1.2]);
%! xlabel ("predictor"); ylabel ("response");
%! title ("actual vs predicted values");
%! legend ({"Predicted Response", "Actual Response", "Theoretical Response"});

## Test output
%!x1 = [-0.87382;-0.66688;-0.54738;0.84894;0.22472;-0.92622;-0.26183;0.55114;-0.7225;-0.24424];
%!x2 = [-0.39428;-0.41286;-0.44929;-0.36923;0.41461;0.46637;-0.63004;0.90951;0.93747;-0.96886];
%!Y  = [-0.92905;-0.48709;-0.16198;-0.87861;0.8525;-0.83343;0.45695;0.66981;0.26195;-0.16609];
%!test
%!ypred = gampredict ([x1, x2], Y, [x1, x2]);
%!yex   = [-0.92889;0.47067;-0.17343;-0.88472;0.85391;-0.83498;0.45797;0.67111;0.26083;-0.16618];
%!assert (ypred, [[-0.92889;0.47067;-0.17343;-0.88472;0.85391;-0.83498;0.45797;0.67111;0.26083;-0.16618]);

## Test input validation
%!error<gampredict: too few input arguments.> gampredict (1)
%!error<gampredict: too few input arguments.> gampredict (1, 2)
%!error<gampredict: number of rows in X and Y must be equal.>  ...
%! gampredict (ones(4,5), ones(5,4), ones(1,1))
%!error<gampredict: number of columns in Xfit must be equal to X.> ...
%! gampredict (ones(4,5), ones(4,1),ones(1,1))
%!error<gampredict: Invalid value in Formula, Formula must be a string.> ....
%! gampredict(ones (5, 5), ones (5, 1), ones (5, 5), 'formula', ones(1, 1))
%!error<gampredict: ResponseName must be a string or char.> ...
%! gampredict(ones (5, 5), ones (5, 1), ones (5, 5), 'responsename', ones(1, 1))
%!error<gampredict: Value in fitstd must be logical or boolean.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'fitstd', 'a')
%!error<gampredict: MaxPValue must be numeric and in between 0 and 1.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'maxpval', 'a')
%!error<gampredict: MaxPValue must be numeric and in between 0 and 1.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'maxpval', 2)
%!error<gampredict: PredictorNames must be supplied as cellstr.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'predictors', 2)
%!error<gampredict: PredictorNames must have same number of columns as X.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'predictors', {'a','b','c'})
%!error<gampredict: Alpha must be numeric and in between 0 and 1.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'Alpha', 'a')
%!error<gampredict: Alpha must be numeric and in between 0 and 1.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'Alpha', 2)
%!error<gampredict: Value in Includeinteractions must be logical or boolean.>...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'includeinteractions', 2)
%!error<gampredict: degree of freedom must be a numeric matrix.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'dof', 2)
%!error<gampredict: degree of freedom must be a numeric matrix.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'dof', [1,2,3])
