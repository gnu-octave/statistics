## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.
##
##
## -*- texinfo -*-
## @deftypefn {statistics} {@var{b} =} ridge (@var{y}, @var{X}, @var{alpha})
## @deftypefnx {statistics} {@var{b} =} ridge (@var{y}, @var{X}, @var{alpha}, @var{flag})
##
## Here,
## Arguments are :
##
## @itemize
## @item
## @code{y} is the predictor matrix or the observed value matrix.
## @item
## @code{X} is the regressor matrix.
## @item
## @code{alpha} is the ridge parameter.
## @item
## @code{flag} is an indicator, pass 0 for restoring the ridge coefficient(s)
## back to the original scale of the data given.
## @end itemize
##
## There are two Types of uses of this function,
##
##  @code{@var{b} = ridge (@var{y}, @var{X}, @var{alpha}) returns the vector of
##  coefficient estimates by applying ridge regression from predictor matrix @var{X}
##  and the response @var{y}. Each value of @var{b} is the coefficient for the
##  respective ridge parameter given @var{alpha}. by default @var{b} is calculated
##  after centering and scaling to mean 0 and standard deviation 1.
##  In general it can be used to produce ridge traces ( see demo example ) with
##  coefficients as the function returns a vector of regression coefficients.
##
##  and,
##
##  @code{@var{b} = ridge (@var{y}, @var{X}, @var{alpha}, @var{flag}) performs the regression
##  and returns the coefficient estimates for specified scaling by @var{flag}. when @var{flag}
##  is given 0,
##  the function restores the coefficient to the scale of the original data thus is more useful for
##  making predictions. for @var{flag} = 1 is same as the above function as no
##  restoration is done to the scale.
##
## Demo,
## @example
## demo ridge
## @end example
## @seealso {lasso, stepwisefit, regress}
## @end deftypefn

function b = ridge (y, X, alpha, flag)

  ## input parameter check
  if  (nargin > 4 || nargin < 3)
    error ("ridge: Invalid number of arguments.");
  endif


  if (! ismatrix (y))
    error ("ridge: Y must be a numeric matrix.");
  endif
  if (! ismatrix (X))
    error ("ridge: X must be a numeric matrix.");
  endif


  if (nargin < 4 || isempty(flag))
    unscale = false;
  elseif ( flag == 1)
    unscale = false;
  elseif ( flag == 0)
    unscale = true;
  endif


  if (columns (y) != 1)
    error ("ridge: Y must be a column vector.");
  endif


  if (rows (y) != rows (X))
    error ("ridge: Y and X must contain the same number of rows.");
  endif

  ## remove any missing values
  notnans = !logical (sum (isnan ([y X]), 2));
  y = y(notnans);
  X = X(notnans,:);

  ## normalising X to mean zero and Std deviation ones
  m = mean(X);
  stdx = std(X,0,1);

  z = (X - m) ./ stdx;

  ## add pseudo observations
  Z_pseudo = [z; (sqrt(alpha(1)) .* eye (columns(X)))];
  Y_pseudo = [y; zeros(columns(X), 1)];
  ## coefficient

  b = Z_pseudo \ Y_pseudo;


  alphas = numel(alpha);


  ## compute the coefficient estimates for additional ridge parameters.
  if (alphas >= 2)

    ## adding a multiple of the identity matrix to the last p rows.

    ## b is set to 0 for the current ridge parameter value
    b(end,alphas) = 0;


   for i=2:alphas
     Z_pseudo ( end - columns(X) + 1 : end,:)  =  sqrt( alpha(i)) .* eye( columns(X));
     b(:,i) = Z_pseudo \ Y_pseudo;
   endfor
  endif

   ## changing back to the scale
   if (unscale)
    b = b ./ repmat (stdx', 1, alphas);
    b = [ mean(y)-m.*b; b];
   endif

endfunction

%!demo
%!
%! load acetylene
%!
%! X = [x1 x2 x3];
%!
%! x0 = ones(16,1);
%! x1x2 = x1 .* x2;
%! x1x3 = x1 .* x3;
%! x2x3 = x2 .* x3;
%!
%! D = [x1 x2 x3 x1x2 x1x3 x2x3];
%!
%! k = 0:1e-5:5e-3;
%!
%! b = ridge(y,D,k);
%!
%! figure
%! plot(k,b,'LineWidth',2)
%! ylim([-100 100])
%! grid on
%! xlabel('Ridge Parameter')
%! ylabel('Standardized Coefficient')
%! title('Ridge Trace')
%! legend('x1','x2','x3','x1x2','x1x3','x2x3')
%!
## test input validation
%!error<ridge: Invalid number of arguments.> ridge(1, 2)
%!error<ridge: function called with too many inputs> ridge(1, 2, 3, 4, 5)
%!error<ridge: Y must be a column vector.> ridge([1, 2, 3], [], 4)
%!error<ridge: Y and X must contain the same number of rows.> ridge ([1; 2; 3; 4; 5], [1, 2, 3], 3)
