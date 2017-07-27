## Copyright (C) 2013-2014 Nir Krakauer <nkrakauer@ccny.cuny.edu>
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
## @deftypefn {Function File} {@var{X_use}, @var{b}, @var{bint}, @var{r}, @var{rint}, @var{stats} =} stepwisefit (@var{y}, @var{X}, @var{penter} = 0.05, @var{premove} = 0.1)
## Linear regression with stepwise variable selection.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{y} is an @var{n} by 1 vector of data to fit.
## @item
## @var{X} is an @var{n} by @var{k} matrix containing the values of @var{k} potential predictors. No constant term should be included (one will always be added to the regression automatically).
## @item
## @var{penter} is the maximum p-value to enter a new variable into the regression (default: 0.05).
## @item
## @var{premove} is the minimum p-value to remove a variable from the regression (default: 0.1).
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{X_use} contains the indices of the predictors included in the final regression model. The predictors are listed in the order they were added, so typically the first ones listed are the most significant.
## @item
## @var{b}, @var{bint}, @var{r}, @var{rint}, @var{stats} are the results of @code{[b, bint, r, rint, stats] = regress(y, [ones(size(y)) X(:, X_use)], penter);}
## @end itemize
## @subheading References
##
## @enumerate
## @item
## N. R. Draper and H. Smith (1966). @cite{Applied Regression Analysis}. Wiley. Chapter 6.
##
## @end enumerate
## @seealso{regress}
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Description: Linear regression with stepwise variable selection

function [X_use, b, bint, r, rint, stats] = stepwisefit(y, X, penter = 0.05, premove = 0.1)

#remove any rows with missing entries
notnans = !any (isnan ([y X]) , 2);
y = y(notnans);
X = X(notnans,:);

n = numel(y); #number of data points
k = size(X, 2); #number of predictors

X_use = [];
v = 0; #number of predictor variables in regression model

iter = 0;
max_iters = 100; #maximum number of interations to do

r = y;
while 1

  iter++;
  #decide which variable to add to regression, if any
  added = false;
  if numel(X_use) < k
    X_inds = zeros(k, 1, "logical"); X_inds(X_use) = 1;
    [~, i_max_corr] = max(abs(corr(X(:, ~X_inds), r))); #try adding the variable with the highest correlation to the residual from current regression
    i_max_corr = (1:k)(~X_inds)(i_max_corr); #index within the original predictor set
    [b_new, bint_new, r_new, rint_new, stats_new] = regress(y, [ones(n, 1) X(:, [X_use i_max_corr])], penter);
    z_new = abs(b_new(end)) / (bint_new(end, 2) - b_new(end));
    if z_new > 1 #accept new variable
      added = true;
      X_use = [X_use i_max_corr];
      b = b_new;
      bint = bint_new;
      r = r_new;
      rint = rint_new;
      stats = stats_new;
      v = v + 1;
    endif
  endif
  
  #decide which variable to drop from regression, if any
  dropped = false;
  if v > 0
    t_ratio = tinv(1 - premove/2, n - v - 1) / tinv(1 - penter/2, n - v - 1); #estimate the ratio between the z score corresponding to premove to that corresponding to penter
    [z_min, i_min] = min(abs(b(2:end)) ./ (bint(2:end, 2) - b(2:end)));
    if z_min < t_ratio #drop a variable
      dropped = true;
      X_use(i_min) = [];
      [b, bint, r, rint, stats] = regress(y, [ones(n, 1) X(:, X_use)], penter);      
      v = v - 1;
    endif
  endif
  
  #terminate if no change in the list of regression variables
  if ~added && ~dropped
    break
  endif

  if iter >= max_iters
    warning('stepwisefit: maximum iteration count exceeded before convergence')
    break
  endif
  
endwhile

endfunction

%!test
%! % Sample data from Draper and Smith (n = 13, k = 4)
%! X = [7 1 11 11 7 11 3 1 2 21 1 11 10; ...
%!     26 29 56 31 52 55 71 31 54 47 40 66 68; ...
%!     6 15 8 8 6 9 17 22 18 4 23 9 8; ...
%!     60 52 20 47 33 22 6 44 22 26 34 12 12]';
%! y = [78.5 74.3 104.3 87.6 95.9 109.2 102.7 72.5 93.1 115.9 83.8 113.3 109.4]';
%! [X_use, b, bint, r, rint, stats] = stepwisefit(y, X);
%! assert(X_use, [4 1])
%! assert(b, regress(y, [ones(size(y)) X(:, X_use)], 0.05))
