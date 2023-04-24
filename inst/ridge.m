## Copyright (C) 2023, Mohammed Azmat alphahan <azmat.dev0@gmail.com>
##
## This file is part of the statistics pacalphaage for GNU Octave.
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
## @deftypefn  {statistics} {@var{b} =} ridge (@var{y}, @var{X}, [@var{alpha}])
## @deftypefnx {statistics} {@var{b} =} ridge (@var{y}, @var{X}, [@var{alpha}], @var{flag})
##
## @code{@var{b} = ridge (@var{y}, @var{X}, @var{alpha}) returns the coefficient
## estimates for ridge ridgeion from ridgeor matrix @var{X} and the response @var{y}.
## Each value of @var{b} is the coefficient for the respective ridge parameter given @var{alpha}.
## by default @var{b} is calculated after centering and scaling to mean 0 and standard deviation 1.
##
## @code{@var{b} = ridge (@var{y}, @var{X}, @var{alpha}, @var{flag}) returns the
## coefficient estimates for specified scaling by @var{flag}. when @var{flag} is given 0,
## the function restores the coefficient to the scale of the original data.
## for @var{flag} = 1 is same as the above function as no restoration is done.
##
## @seealso{ridge, lasso, stepwisefit}
##
## @end deftypefn



function b = ridge (y,X,alpha,flag)

  % input parameter checalpha
  if  (nargin > 4 || nargin < 3)
    error ("ridge: Invalid number of arguments.");
  endif


  if (! ismatrix (y))
    error ("ridge: y must be a numeric matrix.");
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
    error ("ridge: y must be a column vector.");
  endif


  if (rows (y) != rows (X))
    error ("ridge: y and X must contain the same number of rows.");
  endif

  % remove any missing values
  notnans = !logical (sum (isnan ([y X]), 2));
  y = y(notnans);
  X = X(notnans,:);

  % normalising X to mean zero and Std deviation ones
  m = mean(X);
  stdx = std(X,0,1);

  z = (X - m) ./ stdx;

  % add pseudo observations
  Z_pseudo = [X; (sqrt(alpha(1)) * eye (rows(X)))];
  Y_pseudo = [y; zeros(rows(X), 1)];
  % coefficient
  b = Z_pseudo \ Y_pseudo;

  alphas = numel(alpha);

  % compute the coefficient estimates for additional ridge parameters.
  if (alphas >= 2)

    % adding a multiple of the identity matrix to the last p rows.

    % b is set to 0 for the current ridge parameter value
   b(alphas) = 0;

   for i=2:alphas
     Z_pseudo ( end - rows(X) + 1 : end)  =  sqrt( alpha(i) * eye( rows(X)));
     b(:,i) = Z_pseudo \ Y_pseudo;

   endfor
  endif

   % changing back to the scale
   if( unscale)
   b = b ./ repmat( stdx', 1, alphas);
   b = [ mean(y) - m * b; b];
   endif

endfunction


## test input validation
%!error<ridge: Invalid number of arguments.> ridge (1, 2)
%!error<ridge: Invalid number of arguments.> ridge (1, 2, 3, 4, 5)
%!error<ridge: y must be a column vector.> ridge ([1, 2, 3], [], 4)
%!error<ridge: y and X must contain the same number of rows.> ridge ([1; 2; 3; 4; 5], [1, 2, 3], 3)
%!
%!
%!demo
%!
%! load carbig
%! X = [Acceleration Weight Displacement Horsepower];
%! y = MPG;
%!
%! n = length(y);
%! rng('default')
%! c = cvpartition(n,'HoldOut',0.3);
%! idxTrain = training(c,1);
%! idxTest = ~idxTrain;
%!
%! Xt = X(idxTrain,:);
%! yt = y(idxTrain);
%!
%! k = 5;
%!
%! b = ridge(yt,Xt(:,1),k,0);
%!
%! yhat = b(1) + X(idxTest,:)*b(2:end);
%!
%! scatter(y(idxTest),yhat)
%! hold on
%! plot(y(idxTest),y(idxTest))
%! xlabel('Actual MPG')
%! ylabel('Predicted MPG')
%! hold off






