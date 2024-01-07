## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{yFit} =} predict (@var{obj}, @var{Xfit})
## @deftypefnx {statistics} {@var{yFit} =} predict (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{yFit}, @var{ySD}, @var{yInt}] =} predict (@dots{})
##
## Predict new data points using generalized additive model regression object.
##
## @code{@var{yFit} = predict (@var{obj}, @var{Xfit}} returns a vector of
## predicted responses, @var{yFit}, for the predictor data in matrix @var{Xfit}
## based on the Generalized Additive Model in @var{obj}.  @var{Xfit} must have
## the same number of features/variables as the training data in @var{obj}.
##
## @itemize
## @item
## @var{obj} must be a @qcode{RegressionGAM} class object.
## @end itemize
##
## @code{[@var{yFit}, @var{ySD}, @var{yInt}] = predict (@var{obj}, @var{Xfit}}
## also returns the standard deviations, @var{ySD}, and prediction intervals,
## @var{yInt}, of the response variable @var{yFit}, evaluated at each
## observation in the predictor data @var{Xfit}.
##
## @code{@var{yFit} = predict (@dots{}, @var{Name}, @var{Value})} returns the
## aforementioned results with additional properties specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @multitable @columnfractions 0.28 0.02 0.7
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{"alpha"} @tab @tab significance level of the prediction
## intervals @var{yInt}, specified as scalar in range @qcode{[0,1]}. The default
## value is 0.05, which corresponds to 95% prediction intervals.
##
## @item @qcode{"includeinteractions"} @tab @tab a boolean flag to include
## interactions to predict new values based on @var{Xfit}.  By default,
## @qcode{"includeinteractions"} is @qcode{true} when the GAM model in @var{obj}
## contains a @qcode{obj.Formula} or @qcode{obj.Interactions} fields. Otherwise,
## is set to @qcode{false}.  If set to @qcode{true} when no interactions are
## present in the trained model, it will result to an error.  If set to
## @qcode{false} when using a model that includes interactions, the predictions
## will be made on the basic model without any interaction terms.  This way you
## can make predictions from the same GAM model without having to retrain it.
## @end multitable
##
## @seealso{fitrgam, @@RegressionGAM/RegressionGAM}
## @end deftypefn

function [yFit, ySD, yInt] = predict (obj, Xfit, varargin)

  ## Check for sufficient input arguments
  if (nargin < 2)
    error ("@RegressionGAM/predict: too few arguments.");
  endif

  ## Check for valid XC
  if (isempty (Xfit))
    error ("@RegressionGAM/predict: Xfit is empty.");
  elseif (columns (obj.X) != columns (Xfit))
    error (strcat (["@RegressionGAM/predict: Xfit must have the same"], ...
                   [" number of features (columns) as in the GAM model."]));
  endif

  ## Clean Xfit data
  notnansf  = ! logical (sum (isnan (Xfit), 2));
  Xfit      = Xfit (notnansf, :);

  ## Default values for Name-Value Pairs
  alpha = 0.05;
  if (isempty (obj.IntMatrix))
    incInt = false;
  else
    incInt = true;
  endif

  ## Parse optional arguments
  while (numel (varargin) > 0)
    switch (tolower (varargin {1}))

      case "includeinteractions"
        tmpInt = varargin{2};
        if (! islogical (tmpInt) || (tmpInt != 0 && tmpInt != 1))
          error (strcat (["@RegressionGAM/predict: includeinteractions"], ...
                         [" must be a logical value."]));
        endif
        ## Check model for interactions
        if (tmpInt && isempty (obj.IntMat))
          error (strcat (["@RegressionGAM/predict: trained model"], ...
                         [" does not include any interactions."]));
        endif
        incInt = tmpInt;

      case "alpha"
        alpha = varargin{2};
        if (! (isnumeric (alpha) && isscalar (alpha)
                                  && alpha > 0 && alpha < 1))
          error (strcat (["@RegressionGAM/predict: alpha must be a"], ...
                         [" scalar value between 0 and 1."]));
        endif

      otherwise
        error (strcat(["@RegressionGAM/predict: invalid NAME in"], ...
                      [" optional pairs of arguments."]));
    endswitch
    varargin (1:2) = [];
  endwhile

  ## Choose whether interactions must be included
  if (incInt)
    if (! isempty (obj.Interactions))
      ## Append interaction terms to the predictor matrix
      for i = 1:rows (obj.IntMat)
        tindex = logical (obj.IntMat(i,:));
        Xterms = Xfit(:,tindex);
        Xinter = ones (rows (Xfit), 1);
        for c = 1:sum (tindex)
          Xinter = Xinter .* Xterms(:,c);
        endfor
        ## Append interaction terms
        Xfit = [Xfit, Xinter];
      endfor
    else
      ## Add selected predictors and interaction terms
      XN = [];
      for i = 1:rows (obj.IntMat)
        tindex = logical (obj.IntMat(i,:));
        Xterms = Xfit(:,tindex);
        Xinter = ones (rows (Xfit), 1);
        for c = 1:sum (tindex)
          Xinter = Xinter .* Xterms(:,c);
        endfor
        ## Append selected predictors and interaction terms
        XN = [XN, Xinter];
      endfor
      Xfit = XN;
    endif
    ## Get parameters and intercept vectors from model with interactions
    params = obj.ModelwInt.Parameters;
    Interc = obj.ModelwInt.Intercept;
    ## Update length of DoF vector
    DoF = ones (1, columns (Xfit)) * obj.DoF(1);
  else
    ## Get parameters and intercept vectors from base model
    params = obj.BaseModel.Parameters;
    Interc = obj.BaseModel.Intercept;
    ## Get DoF from model
    DoF = obj.DoF;
  endif


  ## Predict values from testing data
  yFit = predict_val (params, Xfit, Interc);

  ## Predict Standard Deviation and Intervals of estimated data (if requested)
  if (nargout > 0)
    ## Ensure that RowsUsed in the model are selected
    Y = obj.Y(logical (obj.RowsUsed));
    X = obj.X(logical (obj.RowsUsed), :);
    ## Predict response from training predictor data with the trained model
    yrs = predict_val (params, X , Interc);

    ## Get the residuals between predicted and actual response data
    rs     = Y - yrs;
    var_rs = var (rs);
    var_pr = var (yFit);  # var is calculated here instead take sqrt(SD)

    t_mul  = tinv (1 - alpha / 2, obj.DoF);

    ydev   = (yFit - mean (yFit)) .^ 2;
    ySD    = sqrt (ydev / (rows (yFit) - 1));

    varargout{1} = ySD;
    if (nargout > 1)
      moe    = t_mul (1) * ySD;
      lower  = (yFit - moe);
      upper  = (yFit + moe);
      yInt   = [lower, upper];

      varargout{2} = yInt;
    endif
  endif

endfunction

## Helper function for making prediction of new data based on GAM model
function ypred = predict_val (params, X, intercept)
  [nsample, ndims_X] = size (X);
  ypred = ones (nsample, 1) * intercept;
  ## Add the remaining terms
  for j = 1:ndims_X
    ypred = ypred + ppval (params(j), X (:,j));
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
%! a = fitrgam (X, Y, "order", [5, 5]);
%! [ypred, ySDsd, yInt] = predict (a, X);
%!
%! ## Plot the results
%! figure
%! [sortedY, indY] = sort (Ytrue);
%! plot (sortedY, "r-");
%! xlim ([0, 80]);
%! hold on
%! plot (ypred(indY), "g+")
%! plot (yInt(indY,1), "k:")
%! plot (yInt(indY,2), "k:")
%! xlabel ("Predictor samples");
%! ylabel ("Response");
%! title ("actual vs predicted values for function f1(x) = cos (3x) ");
%! legend ({"Theoretical Response", "Predicted Response", "Prediction Intervals"});
%!
%! ## Use 30% Holdout partitioning for training and testing data
%! C = cvpartition (80, "HoldOut", 0.3);
%! [ypred, ySDsd, yInt] = predict (a, X(test(C),:));
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


## Test input validation
%!error<@RegressionGAM/predict: too few arguments.> ...
%! predict (RegressionGAM (ones(10,1), ones(10,1)))
%!error<@RegressionGAM/predict: Xfit is empty.> ...
%! predict (RegressionGAM (ones(10,1), ones(10,1)), [])
%!error<@RegressionGAM/predict: Xfit must have the same number of features> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)), 2)
%!error<@RegressionGAM/predict: invalid NAME in optional pairs of arguments.> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)), ones (10,2), "some", "some")
%!error<@RegressionGAM/predict: includeinteractions must be a logical value.> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)), ones (10,2), "includeinteractions", "some")
%!error<@RegressionGAM/predict: includeinteractions must be a logical value.> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)), ones (10,2), "includeinteractions", 5)
%!error<@RegressionGAM/predict: alpha must be a scalar value between 0 and 1.> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)), ones (10,2), "alpha", 5)
%!error<@RegressionGAM/predict: alpha must be a scalar value between 0 and 1.> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)), ones (10,2), "alpha", -1)
%!error<@RegressionGAM/predict: alpha must be a scalar value between 0 and 1.> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)), ones (10,2), "alpha", 'a')
