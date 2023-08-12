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
## @deftypefn  {statistics} {@var{yFit} =} gampredict (@var{X}, @var{y}, @var{Xfit})
## @deftypefnx {statistics} {@var{yFit} =} gampredict (@var{X}, @var{y}, @var{Xfit}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{yFit}, @var{ySD}, @var{yInt}] =} gampredict (@var{X}, @var{y}, @var{Xfit}, @var{name}, @var{value})
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
## @item @tab @qcode{"weights"} @tab Observational weights specified as a
## numeric matrix with each row correspoding to the observations in @qcode{"X"}.
## @qcode{"weights"} must have same number of rows as @qcode{"X"}. Default is
## ones (size (X,1),1).
##
## @item @tab @qcode{"alpha"} @tab Significance level of the prediction
## intervals @qcode{"yInt"}. Specified as scalar in range [0,1]. This argument
## is only valid when @qcode{"fitstd"} is set true. default value is 0.05.
## for example 'alpha',0.05 return 95% prediction intervals.
##
## @item @tab @qcode{"dof"} @tab Degree of freedom to fit a third order spline.
## for fitting a spline @qcode{"dof = knots + order"}, for fitting a GAM a
## polynomial spline of degree '3' is used hence the number of knots can be
## controlled by degree of freedom, degree of freedom can be used to adjust the
## fit of the each variable. the length of @qcode{"dof"} must be same as the
## columns of @qcode{"X"}. default value is 8 for each predictor variable.
##
## @item @tab @qcode{"includeinteractions"} @tab
##
##
## @end multitable
##
## @code{@var{yFit} = knnpredict (@var{X}, @var{y}, @var{Xfit})} returns a
## Numeric Matrix of predicted Response values for predictor data in @var{Xfit},
## using a generalised additive model trained using predictor values in @var{X}
## and response values in @var{y}.
##
##
## @code{@var{yFit} = knnpredict (@var{X}, @var{y}, @var{Xfit}, @var{name}, @var{value})}
## returns a matrix of predicted Response values for predictor data
## in @var{Xfit} using GAM regression model with specified options as Name-value
## pairs.
##
## @code{[@var{yFit}, @var{ySD}, @var{yInt}] = knnpredict (@var{X}, @var{y}, @var{Xfit})}
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
## Each row of the prediction intervals contains @math{100 x (1-alpha) %}
## of prediction interval.
##
##@end multitable
##
## @end itemize
## for demo use demo gampredict
##
## @seealso{regress}
## @end deftypefn

function [yFit, ySD, yInt] = gampredict (X, y, Xfit, varargin)


  ## check positional parameters
  if (nargin < 3)
	  error ("gampredict: too few input arguments.");
  endif

  if (rows (X) != rows (y))
    error ("gampredict: number of rows in X and Y must be equal.");
  endif

  if (columns (X) != columns (Xfit))
    error ("gampredict: number of columns in Xfit must be equal to X.");
  endif

  ## process optional parameters

  ## adding default values
  fitstd  = false;                       ## flag to fit model for std dev
  maxPval = 1;                           ## MAx p-value for detecting interaction terms
  alpha   = 0.05;                        ## significance level for prediction
  catpred = [];                          ## list of categorical predictors
  Interactions   = [];                   ## Nos or list of interaction term
  PredictorNames = [];                   ## predictor variable names
  ResponseName   = "Y";                  ## response variable name
  weights        = ones (size (X,1), 1); ## observational weights
  IncludeInt     = false;                ## for predicting the values
  formulas       = [];                   ## formula for GAM model
  dof            = ones(1,columns(X))*8; ## degree of freedom for fitting spline
  tol            = 1e-3;                 ## positive tolerence to decide convergence


  while (numel (varargin) > 0)
    switch (tolower (varargin {1}))
      case "formula"
        formulas = varargin{2};
      case "responsename"
        ResponseName = varargin{2};
      case "fitstd"
        fitstd = varargin{2};
      case "interactions"
        Interactions = varargin{2};
      case "maxpval"
        maxPval = varargin{2};
      case "catpredictors"
        catpred = varargin{2};
      case "predictors"
        PredictorNames = varargin{2};
      case "alpha"
        alpha = varargin{2};
      case "weights"
        weights = varargin{2};
      case "includeinteractions"
        IncludeInt = varargin{2};
      case "dof"
        dof = varargin{2};
      otherwise
        error ("gampredict: invalid NAME in optional pairs of arguments.")
    endswitch
    varargin (1:2) = [];
  endwhile

  ##------checking optional parameters------##

  ## formulas
  if (! isempty (formulas))
    if (! ischar (formulas))
      error ("gampredict: Invalid value in Formula, Formula must be a string.");
    endif
  endif

  ## ResponseName
  if (! ischar (ResponseName))
    error ("gampredict: ResponseName must be a string or char.");
  endif

  ## fitstd
  if (! islogical (fitstd) && fitstd != 1 && fitstd != 0)
    error ("gampredict: Value in fitstd must be logical or boolean.");
  endif

  ## Interactions
  if (! isempty (Interactions))
    if (! islogical (Interactions) || ! ismatrix (Interactions))
      error ("gampredict: Invalid Value(s) in Interactions.");
    endif
  endif


  ## maxPval
  if (! isscalar (maxPval) || maxPval < 0 || maxPval > 1)
    error ("gampredict: maxPval must be numeric and in between 0 and 1.");
  endif

  ## catpred

  ## predictorNames
  if (! isempty (PredictorNames))
    if (! iscellstr (PredictorNames))
      error ("gampredict: PredictorNames must be supplied as cellstr.");
    elseif (columns (PredictorNames) != columns (X))
      error ("gampredict: PredictorNames must have same number of columns as X.");
    endif
  else
    ## empty predictornames generate default
    PredictorNames = gendefault (columns (X));
  endif

  ## Alpha
  if (! isscalar (alpha) || alpha < 0 || alpha > 1)
    error ("gampredict: alpha must be numeric and in between 0 and 1.");
  endif

  ##includeinteractions
  if (! islogical (IncludeInt) && IncludeInt != 1 && IncludeInt != 0)
    error ("gampredict: Value in Includeinteractions must be logical or boolean.");
  endif

  ## DF
  if (! ismatrix (dof) || ! ismatrix (dof) || columns (dof) != columns (X))
    error ("gampredict: degree of freedom must be a numeric matrix.");
  endif

  ##------optional parameter check end -----##

  ## process predictor matrix
  notnans  = ! logical (sum (isnan ([y, X]), 2));
  notnansf = ! logical (sum (isnan (Xfit), 2));
  y        = y (notnans);
  X        = X (notnans, :);
  Xfit     = Xfit (notnansf, :);



  ## fit the GAM model

  order     = 3;
  knots     = dof - order;  ## from definition
  intercept = mean (y);
  res       = y - intercept;
  [ppfk, RSS] = onecycleBackfit (X, y, res, knots (1), intercept, tol);
  converged = false;


  ## main loop
  do
    ## single cycle of backfit
    for j = 1:columns (X)

      ## calculate residuals to fit spline
      res = res + intercept + ppval (ppfk (j), X (:, j));
      gk = smoother (X (:,j), res, knots (j), order);

      ## centering coeffs
      gk.coefs = gk.coefs - sum (sum (gk.coefs)) / rows (X);
      RSSk (j) = abs (sum (abs (y - ppval (gk, X(:,j)) - intercept )) .^ 2) / rows (y);
      ppfk (j) = gk;
      res = res - intercept - ppval (ppfk (j), X (:,j));
    endfor
    RSS - RSSk

    ## check if RSS is less than the tolerence
    if (all (abs (RSS - RSSk) <= tol))
      converged = true;
    endif

    ## update RSS
    RSS = RSSk;
  until (converged)


  ## predict values
  yFit = predict_val (ppfk, Xfit, intercept);

  ## calculate yInt
  if (fitstd)
    rs     = yFit - mean(y);
    var_rs = var (rs);
    var_pr = var (yFit);
    t_mul  = tinv (1 - alpha / 2, dof);
    ytru   = sort (yFit);
    moe    = t_mul * sqrt (var_pr + var_rs);

    lower  = (yFit - moe);
    upper  = (yFit + moe);

    yInt   = [lower, upper];


  endif

endfunction

function pp = smoother (x, res, knots, order)
  ## sort the data and pass it to splinefit to fit a spline
  pp = splinefit (x, res, knots, 'order', order-1);
endfunction

## ------ predict_val -------- ##
function ypred = predict_val (splstr, X, intercept)

  [n_samples n_features] = size (X);
  ypred = ones (rows (X), 1) * intercept;

  for j = 1:n_features
    ypred = ypred + ppval (splstr(j), X (:,j));
  endfor
endfunction
## ------ predict_val end ------ ##

## Function parseInteractions to detect user given interaction
function intMat = parseInteractions (formulas, predictorNames)
  intMat = [];
  if (islogical (formulas))
    if (numel (predictorNames) != columns (formulas))
      error ("gampredict: ");
    endif
    for i = 1:rows (formulas)
      ind = find (formulas(i, :) == 1);
      intMat = [intMat; ind];
    endfor
  elseif (ischar (formulas))
    if (strcmpi (formulas, "all"))
      #calculate all p*(p-1)/2 interaction terms

    else
      #calculate interaction matrix by formula
      formulaParts = strsplit(formulas, '~');
      responseVar = strtrim(formulaParts{1});
      predictorString = strtrim(formulaParts{2});
      intterms = strsplit(predictorString, '+');
      in = [];
      for i = 1:numel (intterms)
        if (cell2mat (strfind (intterms (i), ':')))
          inters = strtrim (strsplit (cell2mat (intterms (i)), ':'));
          in = [in; (index (predictorNames, inters(1)) + ...
                                          index (predictorNames, inters(2)))];
        endif
      endfor
      ## pass recursively to get intMat
      intMat = parseInteractions (logical(in), predictorNames);
    endif
  else
    error ("gampredict: Invalid value in Interactions.");
  endif
endfunction
##------ parseInteractions End ------##

##------ gendefault to generate default predictornames-----##
function defs = gendefault (p)
  for i = 1:p
    defs {i} = strcat ('x', num2str (i));
  endfor
endfunction
## ------gendefault end ------ ##

## ------once cycle backfit ------##
function [ppf, RSSk] = onecycleBackfit (X, y, res, knots, m, tol)
for i = 1:2
  for j = 1:columns (X)
    if (i > 1)
      res = res + m + ppval (ppf (j), X (:, j));
    endif

    gk = splinefit (X (:,j), res, knots, 'order', 2);
    gk.coefs = gk.coefs - sum (sum (gk.coefs)) / rows (X);
    RSSk = (sum (abs (y - ppval (gk, X(:,j)) - m )) .^ 2) / rows (y);
    ppf (j) = gk;
    res = res - m - ppval (ppf (j), X (:,j));
    #ppf (j).coefs = ppf (j).coefs - sum (sum (ppf (j).coefs)) ./ rows (X);
  endfor
endfor
endfunction
## ------ onecycleBackfit end -----##



%!demo
%! ## synthetic datasample consisting of two different function
%! f1 = @(x) cos (3 *x);
%! f2 = @(x) x .^ 3;
%!
%! samples = 100;
%! rand("seed",9);
%!
%! ## generating random samples for f1 and f2
%! x1 = 2 * rand (samples, 1) - 1;
%! x2 = 2 * rand (samples, 1) - 1;
%!
%! y = f1(x1) + f2(x2);
%!
%! ## adding noise to the data
%! y = y + y .* 0.2 .* rand (samples,1);
%!
%! X = [x1, x2];
%!
%! ypred = gampredict (X, y, X);
%!
%!
%! subplot (2, 2, 2, "align")
%! plot (x1, ypred, 'o', 'color', 'r', x1, y, 'o', 'color', 'b')
%! xlabel ("x1"); ylabel ("Y");
%! title ("actual vs predicted values for function f1(x) = cos (3x) ");
%! legend ({"predicted by GAM", "Actual Value"});
%!
%! subplot (2, 2, 4, "align")
%! plot (x2, ypred, 'o', 'color', 'r', x2, y, 'o', 'color', 'b')
%! xlabel ("x1"); ylabel ("Y");
%! title ("actual vs predicted values for function f1(x) = x^3 ");
%! legend ({"predicted by GAM", "Actual Value"});
%!
%! subplot (1, 2, 1, "align")
%! plot (y, ypred, 'o')
%! xlabel ("y"); ylabel ("Y_predicted");
%! title ("actual vs predicted values");
%! legend ({"predicted by GAM", "Actual Value"});

## Test output
%!x1 = [-0.87382;-0.66688;-0.54738;0.84894;0.22472;-0.92622;-0.26183;	...
%!     -0.55114;-0.7225;-0.24424;-0.24996;-0.97293;0.90399;0.43813;	 ...
%!     -0.12081;-0.15285;-0.44294;0.16028;-0.44301;0.12117];
%!x2 = [-0.39428;-0.41286;-0.44929;-0.36923;0.41461;0.46637;-0.63004;	...
%!      0.90951;	0.93747;	-0.96886;	0.023887;-0.12016;0.64339;-0.13052;	...
%!      0.15532;	-0.14205;	0.66689;	0.25162;	0.82184;	0.014455];
%!
%!y  = [-0.92905;	-0.48709;	-0.16198;	-0.87861;	0.8525;	-0.83343;	0.45695;	...
%!      0.66981;	 0.26195;	-0.16609;	0.73179;	-0.97701;	-0.6428;	0.2514;	 ...
%!      0.93878;	0.89383;	0.5362;	0.90254;	0.79453;	0.93465];
%!test
%!ypred = gampredict ([x1, x2], y, [x1, x2]);
%!yex   = [-1.3466;-0.79371;-2.3899;1.1113;-2.0822;0.15497;0.83367; ...
%!        0.0036363;-0.72754;0.99625;-2.4546;-1.9636;0.026161;1.4654; ...
%!        1.1415;0.73443;1.272;1.054;1.4404];
%!assert (ypred, yex)

## Test input validation
%!error<gampredict: too few input arguments.> gampredict (1)
%!error<gampredict: too few input arguments.> gampredict (1, 2)
%!error<gampredict: number of rows in X and Y must be equal.>  ...
%! gampredict (ones(4,5), ones(5,4),ones(1,1))
%!error<gampredict: number of columns in Xfit must be equal to X.> ...
%! gampredict (ones(4,5), ones(4,1),ones(1,1))
%!error<gampredict: Invalid value in Formula, Formula must be a string.> ....
%! gampredict(ones (5, 5), ones (5, 1), ones (5, 5), 'formula', ones(1, 1))
%!error<gampredict: ResponseName must be a string or char.> ...
%! gampredict(ones (5, 5), ones (5, 1), ones (5, 5), 'responsename', ones(1, 1))
%!error<gampredict: Value in fitstd must be logical or boolean.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'fitstd', 'a')
%!error<gampredict: maxPval must be numeric and in between 0 and 1.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'maxpval', 'a')
%!error<gampredict: maxPval must be numeric and in between 0 and 1.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'maxpval', 2)
%!error<gampredict: PredictorNames must be supplied as cellstr.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'predictors', 2)
%!error<gampredict: PredictorNames must have same number of columns as X.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'predictors', {'a','b','c'})
%!error<gampredict: alpha must be numeric and in between 0 and 1.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'alpha', 'a')
%!error<gampredict: alpha must be numeric and in between 0 and 1.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'alpha', 2)
%!error<gampredict: Value in Includeinteractions must be logical or boolean.>...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'includeinteractions', 2)
%!error<gampredict: degree of freedom must be a numeric matrix.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'dof', 2)
%!error<gampredict: degree of freedom must be a numeric matrix.> ...
%! gampredict(ones (5,5), ones (5,1), ones (5,5), 'dof', [1,2,3])
