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
## @deftypefn  {statistics} {@var{obj} =} fitrgam
## @deftypefnx {statistics} {@var{obj} =} fitrgam (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} fitrgam (@var{X}, @var{Y}, @var{name}, @var{value})
##
## Create a Generalised additive model (GAM) Regression Model,
## RegressionGAM Object using @var{X}, @var{Y} and other additional 
## Name-Value pairs. Returned Object can be used to predict new
## values, Properties of the Object can be altered via returned Object 
## @var{obj}.
## 
## Input arguments that can be given to create and object of 
## class @qcode{RegressionGAM} are :
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
##
## @emph{Additional parameters can be passed as name value pairs :}
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @tab @qcode{"Xfit"} @tab  must be a @math{MxP} numeric matrix of 
## query/new points to predict the response value.
## @var{Xfit} must have same numbers of columns as @var{X}.
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
## @end multitable
## 
##
##
## @end itemize
## for demo use demo RegressionGAM
##
## @seealso{regress}
## @end deftypefn

function obj = fitrgam (X, Y, varargin)
  
  ## check the nargins if within range
  ## sanity check will be done by constructor
  
  if ( nargin < 2 && nargin != 0)
    error ("fitrgam: Too few arguments.");
  endif
  
  if (nargin > 12)
    error ("fitrgam: Too many arguments.");
  endif
  
  if (nargin == 0)
    ## return an empty object with warning
    obj = RegressionGAM ();
    warning ("fitrgam: No arguments Provided, Created object will be Empty.");
    
  else
    ## arguments within range and not empty
    obj = RegressionGAM (X, Y, varargin {:});
    
  endif
endfunction

%!demo
%! # Train a RegressionGAM Model for synthetic values
%!
%! f1 = @(x) cos (3 *x);
%! f2 = @(x) x .^ 3;
%!
%! # generate x1 and x2 for f1 and f2
%! x1 = 2 * rand (50, 1) - 1;
%! x2 = 2 * rand (50, 1) - 1;
%!
%! # calculate y
%! y = f1(x1) + f2(x2);
%!
%! # add noise
%! y = y + y .* 0.2 .* rand (50,1);
%! X = [x1, x2];
%!
%! # create an object
%! a = fitrgam (X, y, "tol", 1e-3)


## Test input validation
%!error<fitrgam: Too few arguments.> fitrgam (ones(10,2))
%!error<fitrgam: Too many arguments.> ...
%! fitrgam (1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9)
%!warning<fitrgam: No arguments Provided, Created object will be Empty.> ...
%! fitrgam ();
