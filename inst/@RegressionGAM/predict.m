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
##
## -*- texinfo -*-
## @deftypefn  {statistics} {@var{yFit} =} predict (@var{obj})
## @deftypefn  {statistics} {@var{yFit} =} predict (@var{obj}, @var{Xfit})
## @deftypefn  {statistics} {@var{yFit} =} predict (@var{obj}, @var{Xfit}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{yFit}, @var{ySD}, @var{yInt}] =} predict (dots())
##
## Predict new data points using trained GAM model object of class 
## RegressionGAM, by GAM regression
##
## @itemize
## @item
## @code{obj} must be an object of class @code{ClassificationKNN}.
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @tab @qcode{"Alpha"} @tab Significance level of the prediction
## intervals @qcode{"yInt"}. Specified as scalar in range [0,1]. This argument
## is only valid when @qcode{"fitstd"} is set true. default value is 0.05.
## for example 'alpha',0.05 return 95% prediction intervals.
##
## @item @tab @qcode{"includeinteractions"} @tab flag to indicate weather to 
## include interactions to predict new values from.
## @end multitable
## @end itemize
## @end deftypefn
##

function [yFit, ySD, yInt] = predict (obj, Xfit, varargin)
  
  if (nargin < 2 || nargin > 4)
    error ("@RegressionGAM.predict: Too few arguments.");
  endif
  
  ##check obj
  if (! strcmpi (class (obj), "RegressionGAM"))
    error (strcat (["@RegressionGAM.predict: obj must be of"], ...
                   [" Class RegressionGAM."]));
  endif  
  if ( isempty (Xfit))
    error ("@RegressionGAM.predict: Xfit is empty.");
  elseif (columns (obj.X) != columns (Xfit))
    error (strcat ( ["@RegressionGAM.predict: number of columns in Xfit"], ...
                    [" must be equal to X in object."]));
  else
    ## clean Xfit data
    notnansf  = ! logical (sum (isnan (Xfit), 2));
    Xfit      = Xfit (notnansf, :);
  endif
  
  ## Default value for Name-Value Pairs
  Alpha = 0.05;
  includeint = false;
  
  while (numel (varargin) > 0)
    switch (tolower (varargin {1}))
      
      case "includeinteractions"
        includeint = varargin{2};
        if (! islogical (includeint) || includeint != 0 && includeint != 0)
          error (strcat (["@RegressionGAM.predict: includeinteractions"], ...
                         [" must be a logical value."]));
        endif
          
      case "alpha"
        Alpha = varargin{2};
        if (! (isnumeric (Alpha) && isscalar (Alpha) 
                                 && Alpha > 0 && Alpha < 1))
          error (strcat (["@RegressionGAM.predict: alpha must be numeric"], ...
                         [" and in between 0 and 1."]));
        endif
                                     
      otherwise
        error (strcat(["@RegressionGAM.predict: invalid NAME in"], ...
                      [" optional pairs of arguments."]));
    endswitch
    varargin (1:2) = [];
  endwhile
   
  ## Predict values from testing data
  yFit = predict_val (obj.ppfk, Xfit, obj.Intercept);
  
  ## Predict Standard Deviation and Intervals of estimated data (if requested)
  if (nargout > 0)
    ## Predict response from training predictor data with the trained model
    yrs = predict_val (obj.ppfk, obj.X , obj.Intercept);
    ## Get the residuals between predicted and actual response data

    rs     = obj.Y - yrs;
    var_rs = var (rs);    ##
    var_pr = var (yFit);  ## var is calculated here instead take sqrt(SD)

    t_mul  = tinv (1 - Alpha / 2, obj.DoF);


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
%! a = RegressionGAM (X, Y, "order", [5, 5]);
%! [ypred, ySDsd, yInt] = predict (a, X);
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

## Test output

## Test input validation
%!error<@RegressionGAM.predict: Too few arguments.> predict()
%!error<@RegressionGAM.predict: Too few arguments.> predict(1)
%!error<@RegressionGAM.predict: Too few arguments.> predict(1,2,3,4,5)
%!error<@RegressionGAM.predict: obj must be of Class RegressionGAM.> ...
%! predict (1,2)
%!error<@RegressionGAM.predict: Xfit is empty.> ...
%! predict (RegressionGAM (ones(10,1), ones(10,1)),[])
%!error<@RegressionGAM.predict: number of columns in Xfit must be equal to X in object.> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)),2)
%!error<@RegressionGAM.predict: invalid NAME in optional pairs of arguments.> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)), ones (10,2), "some", "some")
%!error<@RegressionGAM.predict: includeinteractions must be a logical value.> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)), ones (10,2), "includeinteractions", "some")
%!error<@RegressionGAM.predict: includeinteractions must be a logical value.> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)), ones (10,2), "includeinteractions", 5)
%!error<@RegressionGAM.predict: alpha must be numeric and in between 0 and 1.> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)), ones (10,2), "alpha", 5)
%!error<@RegressionGAM.predict: alpha must be numeric and in between 0 and 1.> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)), ones (10,2), "alpha", -1)
%!error<@RegressionGAM.predict: alpha must be numeric and in between 0 and 1.> ...
%! predict (RegressionGAM(ones(10,2), ones(10,1)), ones (10,2), "alpha", 'a')
