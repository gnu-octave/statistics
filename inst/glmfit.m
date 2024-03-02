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
## @deftypefn  {statistics} {@var{b} =} glmfit (@var{X}, @var{y}, @var{distribution})
## @deftypefnx {statistics} {@var{b} =} glmfit (@var{X}, @var{y}, @var{distribution}, 
## @var{Name}, @var{Value})
##
## Perform generalized linear model fitting.
##
## This function fits a generalized linear model (GLM) to the given data using
## the specified link function and distribution of the response variable. The
## model is fitted using Iteratively Reweighted Least Squares (IRLS).
##
## @itemize
## @item @var{X} is an @math{nxp} matrix of predictor variables with
## @math{n} observations and @math{p} predictors.
## @item @var{y} is an @math{nx1} vector of response variables.
## @item @var{distribution} specifies the distribution of the response variable
## (e.g., 'poisson').
## @end itemize
##
## @code{@var{b} = glmfit (@var{X}, @var{y}, @var{distribution}, @var{Name},
## @var{Value})} specifies additional options using @qcode{Name-Value} pair
## arguments.
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @qcode{"link"} @tab @tab A character vector specifying a link
## function.
## @end multitable
##
## The function returns @var{b}, the estimated coefficients of the GLM,
## including the intercept term as the first element of @var{b}.
##
## Currently, the function supports only the 'poisson' distribution
## and 'log' link function. Further expansion is required to support
## additional distributions and link functions.
##
## @end deftypefn

function b = glmfit(X, y, distribution, varargin)
  ## Check input
  y = round (y(:));
  if (nargin < 2)
    X = zeros (length (y), 0);
  endif
  xymissing = (isnan (y) | any (isnan (X), 2));
  y(xymissing) = [];
  X(xymissing,:) = [];
  [my, ny] = size (y);
  [mx, nx] = size (X);
  if (mx != my)
    error ("glm: X and y must have the same number of observations.");
  endif

  ## Add defaults
  link = "log";

  ## Parse extra parameters
  while (numel (varargin) > 0)
    switch (tolower (varargin {1}))
      case "link"
        link = varargin{2};
        if (! any (strcmpi (link, {"log"})))
          error ("glmfit: unsupported link function.");
        endif
      otherwise
        error ("glmfit: unknown parameter name.");
    endswitch
    varargin (1:2) = [];
  endwhile

  ## Add column of ones
  X = [ones(size(X, 1), 1), X];
  ## Initialize b
  b = zeros (size (X, 2), 1);
  max_itr  = 1000;
  tolerance = 1e-6;
  b_prev = b;

  ## Select functions
  switch (link)
    case "log"
      inverse_link_func = @(linear_predictor) exp(linear_predictor);
  endswitch
  switch (distribution)
    case "poisson"
      working_response_func = ...
        @(y, linear_predictor, mu) linear_predictor + (y - mu) ./ mu;
      diag_matrix_func = @(mu) diag(mu);
  endswitch


  for i = 1:max_itr
    linear_predictor = X * b;
    ## Give inverse function according to link function
    mu = inverse_link_func (linear_predictor);
    ## Weights for IRLS
    z = working_response_func (y, linear_predictor, mu);
    W = diag_matrix_func (mu);
    ## Update b
    b = (X' * W * X) \ (X' * W * z);
    ## Check for convergence
    if norm(b - b_prev, 2) < tolerance
      break;
    endif
    b_prev = b;
  endfor

  if (i == max_itr)
    #warning('glmfit: reached limit');
  endif
endfunction

## Test output
%!test
%! rand ("seed", 1);
%! X = rand (50, 1);
%! b_true = [0.4; 1.5];
%! mu_true = exp (b_true(1) + b_true(2) * X);
%! y = poissrnd (mu_true);
%! b = glmfit(X, y, "poisson", "link", "log");
%! assert(b(1), b_true(1), 0.5);
%! assert(b(2), b_true(2), 0.5);

## Test input validation
%!error glmfit()
%!error glmfit(rand(5,2))
%!error glmfit(rand(5,2),rand(5,1))
%!error glmfit(rand(5,2), rand(5,1), 'poisson', 'link')
%!error <X must be a numeric.> glmfit('abc', rand(6,1), 'poisson')
%!error <y must be a numeric.> glmfit(rand(5,2), 'abc', 'poisson')
%!error <invalid distribution.> glmfit(rand(5,2), rand(5,1), 2)
%!error <invalid parameter name.> glmfit(rand(5,2), rand(5,1), 'poisson', 2, 'log')
%!error <invalid link function.> glmfit(rand(5,2), rand(5,1), 'poisson', 'link', 2)
%!error <glmfit: X and y must have same number of observations.> glmfit(rand(5,2), rand(6,1), 'poisson')
%!error <glmfit: y cannot be a matrix.> glmfit(rand(10,2), ones(2,2), 'poisson')
%!error <glmfit: unsupported link function.> glmfit(rand(10,2), rand(10,1), 'poisson', 'link', 'inverse')
%!error <glmfit: unsupported parameter name.> glmfit(rand(10,2), rand(10,1), 'poisson', 'notALink', 'log')
%!error <glmfit: distribution must be a recognized value.> glmfit(rand(10,2), rand(10,1), 'xyz')
