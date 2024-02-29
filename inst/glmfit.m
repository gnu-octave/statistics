## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
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

## -*- texinfo -*-
## @deftypefn {statistics} {@var{beta} =} glmfit (@var{X}, @var{y}, @var{distribution}, @var{link})
## Perform generalized linear model fitting 
## 
## This function fits a generalized linear model (GLM) to the given data using 
## the specified link function and distribution of the response variable.
## The model is fitted using Iteratively Reweighted Least Sqaures (IRLS).
##
## @var{X} is an @var{n}-by-@var{p} matrix of predictor variables with @var{n} observations and @var{p} predictions. 
## @var{y} is an @var{n}-by-1 vector of response varibales.
## @var{distribution} specifies the distribution of the response variable (e.g., 'poisson'). 
## @var{link} specifies the link function to use (e.g., 'log').
##
## The function returns @var{beta}, the estimated coefficients of the GLM, including 
## the intercept term as the first element of @var{beta}.
##
## Currently, the function supports only the 'poisson' distribution and 'log' link function.
## Further expansion is required to support additional distributions and link functions.
## 
## @end deftypefn

function beta = glmfit(X, y, distribution, link)
    ## check input
    y = round (y(:));
    if (nargin < 2)
        X = zeros (length (y), 0);
    endif;
    xymissing = (isnan (y) | any (isnan (X), 2));
    y(xymissing) = [];
    X(xymissing,:) = [];
    [my, ny] = size (y);
    [mx, nx] = size (X);
    if (mx != my)
        error ("glm: X and y must have the same number of observations");
    endif

    ## add column of ones
    X    = [ones(size(X, 1), 1), X]; 
    ## initialize beta
    beta = zeros(size(X, 2), 1);  
    max_itr  = 1000;
    tolerance = 1e-6;
    beta_prev = beta;
  
    for i = 1:max_itr
        linear_predictor = X * beta;
        ## give inverse function according to link function
        mu = inverse_link_function(linear_predictor, link);
        ## weights for IRLS 
        z = working_response(distribution, y, linear_predictor, mu);
        W = diag_matrix(distribution, mu);
        ## Update beta 
        beta = (X' * W * X) \ (X' * W * z);
        ## check for convergence
        if norm(beta - beta_prev, 2) < tolerance
            break;
        endif
        beta_prev = beta;
    endfor

    if (i == max_itr)
        disp('warning: reached limit');
    endif
end

function mu = inverse_link_function(linear_predictor, link)
    switch link
        case 'log'
            mu = exp(linear_predictor);
        ## Add cases for other link functions 
        otherwise
            error('Unsupported link function');
    endswitch
end

function z = working_response(distribution, y, linear_predictor, mu)
    switch distribution
        case 'poisson'
            z = linear_predictor + (y - mu) ./ mu;
        ## add other cases
        otherwise
            error('Unsupported distribution');
    endswitch
end

function W = diag_matrix(distribution, mu)
    switch distribution
        case 'poisson'
            w = mu;
            W = diag(w);
        ## add other cases
        otherwise
            error('Unsupported distribution');
    endswitch
end

%!test
%! N = 50;
%! X = rand(N, 1);
%! beta_true = [0.4; 1.5]; 
%! mu_true = exp(beta_true(1) + beta_true(2) * X);
%! y = poissrnd(mu_true);
%! distribution = 'poisson';
%! link = 'log';
%! beta = glmfit([X], y, distribution, link);
%! assert(beta(1), beta_true(1), 0.1);
%! assert(beta(2), beta_true(2), 0.1);

%!test
%! X = (1:20)';
%! beta_true = [1; 0.1];
%! mu_true = exp(beta_true(1) + beta_true(2) * X);
%! y = exp(mu_true);
%! distribution = 'poisson';
%! link = 'log';
%! beta = glmfit(X, y, distribution, link);
%! assert(beta(1), beta_true(1), 0.1);
%! assert(beta(2), beta_true(2), 0.1);