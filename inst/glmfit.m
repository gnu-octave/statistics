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
## @deftypefnx  {statistics} {[@var{b}, @var{dev}] =} glmfit (@var{X}, @var{y}, 
## @var{distribution})
## @deftypefnx  {statistics} {[@var{b}, @var{dev}] =} glmfit (@var{X}, @var{y}, 
## @var{distribution}, 
## @var{Name}, @var{Value})
##
## Perform generalized linear model fitting.
##
## This function fits a generalized linear model (GLM) to the given data using
## the specified link function and distribution of the response variable. 
##
## @itemize
## @item @var{X} is an @math{nxp} matrix of predictor variables with
## @math{n} observations and @math{p} predictors.
## @item @var{y} can be:
##   @itemize
##   @item An @math{n x 1} vector of responses for all the distributions.
##   @item An @math{n x 2} matrix where the first column contains the number of successes
##    and the second column contains the number of trials for the 'binomial' distribution.
##   @end itemize
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
## You can also specify a custom link as a structure with three fields:
## Link Function, Derivative Function, Inverse Function   
## @item @qcode{"constant"} @tab @tab Specifies whether to 
## include a constant term in the model. Options are 
## @qcode{"on"} (default) or @qcode{"off"}.
## @end multitable
## 
## The function returns @var{b}, the estimated coefficients of the GLM,
## including the intercept term as the first element of @var{b}.
##
## @code{@var{[b, dev]} = glmfit (@var{X}, @var{y}, @var{distribution}, @dots{})} 
## also returns the deviance of the fit.
##
## Supported distributions include 'poisson', 'binomial', and 'normal'. 
## The function also supports custom link functions provided as function handles.
##
## Examples:
##
## @example
## @group
## X = [1, 2, 3; 1, 3, 5; 1, 5, 7];
## y = [1; 2; 3];
## distribution = 'poisson';
## [b, dev] = glmfit(X, y, distribution, 'link', 'log');
## @end group
## @end example
##
## @end deftypefn

function [b,varargout] = GLM(X, y, distribution, varargin)
  ## Check input
  if nargin < 3
    error ("glmfit: at least three input arguments (X, y, distribution) are required.");
  endif
  xymissing = (isnan (y) | any (isnan (X), 2));
  y(xymissing) = [];
  X(xymissing,:) = [];
  [my, ny] = size (y);
  [mx, nx] = size (X);
  if (mx != my)
    error ("glmfit: X and y must have the same number of observations.");
  endif

  ## Check dimensions based on distribution
  if strcmpi(distribution, 'binomial')
    if size(y, 2) > 2
      error("glmfit: for a binomial distribution, y must be an n-by-1 or n-by-2 matrix.");
    endif
  else
    if size(y, 2) != 1
      error("glmfit: for non-binomial distributions, y must be an n-by-1 vector.");
    endif
  endif

  ## Add defaults
  switch (tolower (distribution))
    case "poisson"
      link = "log";
    case "binomial"
      link = "logit";
    case "normal"
      link = "identity";
  endswitch

  ## Default
  constant = "on";

  ## Parse extra parameters
  while (numel (varargin) > 0)
    switch (tolower (varargin {1}))
      case "link"
        linkInput = varargin {2};
        ## Check custom link 
        if iscell (linkInput) && numel (linkInput) == 3
          linkFunc = linkInput{1};
          derLinkFunc = linkInput{2};
          invLinkFunc = linkInput{3};
          ## Check for function_handle
          if ~all (cellfun (@(f) isa (f, 'function_handle'), linkInput))
            error ("glmfit: custom link functions must be function handles.")
          endif
          ## Test the custom functions with a small vector
          try
            testInput = rand(5, 1);
            testOutput = invLinkFunc(testInput);
            if ~isequal (size (testInput), size (testOutput))
              error("glmfit: custom inverse link function must return output of the same size as input.");
            endif
          catch
            error("glmfit: error testing custom inverse link function.");
          end_try_catch
          link = "custom";
        ## Check link 
        elseif ischar (linkInput) || isstring (linkInput)
          link = tolower (linkInput);
          if (! any (strcmpi (link, {"identity", "log", "logit", "probit", "loglog", "comploglog", "reciprocal"})))
            error ("glmfit: unsupported link function.");
          endif
        else
          error ("glmfit: link specification is incorrect.");
        endif
      case "constant"
        constant = tolower (varargin {2});
        if (! any (strcmpi (constant, {"on", "off"})))
          error ("glmfit: constant should be set 'on' or 'off'.");
        endif
      otherwise
        error ("glmfit: unknown parameter name.");
    endswitch
    varargin (1:2) = [];
  endwhile

  ## Adjust X based on constant
  if (any (strcmpi (constant, {"on"})))
    X = [ones(size (X, 1), 1), X];
  endif

  ## Initialize b
  b = zeros (size (X, 2), 1);

  ## Select functions
  switch (link)
    case "identity"
      inverse_link_func = @(linear_predictor) linear_predictor;
    case "log"
      inverse_link_func = @(linear_predictor) exp (linear_predictor);
    case "logit"
      inverse_link_func = @(linear_predictor) exp (linear_predictor) ./ (1 + exp (linear_predictor));
    case "probit"
      inverse_link_func = @(linear_predictor) normcdf (linear_predictor);
    case "loglog"
      inverse_link_func = @(linear_predictor) exp (-exp (linear_predictor));
    case "comploglog"
      inverse_link_func = @(linear_predictor) 1 - exp (-exp (linear_predictor));
    case "reciprocal"
      inverse_link_func = @(linear_predictor) 1 ./ linear_predictor;
    case "custom"
      inverse_link_func = invLinkFunc;
  endswitch

  ## Select negLogLikelihood according to distribution
  switch (tolower (distribution))
    case "poisson"
      negLogLikelihood = ...
        @(b) - sum ((y .* X * b) - inverse_link_func (X * b) - gammaln (y+1));
    case "binomial"
      eps = 1e-6;
      if size (y, 2) == 1
        successes = y;
        trials = ones(size (y));
      elseif size (y, 2) == 2
        successes = y(:, 1);
        trials = y(:, 2);
      else
        error("glmfit: For 'binomial' distribution, y must be an n-by-1 or n-by-2 matrix.");
      endif
      negLogLikelihood = @(b) ...
          -sum (successes .* log (max (min (inverse_link_func (X * b), 1 - eps), eps)) ...
          + (trials - successes) .* log (1 - max (min (inverse_link_func (X * b), 1 - eps), eps)));
    case "normal"
      negLogLikelihood = @(b) 0.5 * sum ((y - inverse_link_func (X * b)) .^ 2);
    otherwise
      error("glmfit: unsupported distribution.")
  endswitch

  options = optimset ('MaxFunEvals', 10000, 'MaxIter', 10000);
  b = fminsearch(negLogLikelihood, b, options);
  
  if (nargout > 1)
      dev = [];
      switch (tolower (distribution))
        case "poisson"
          p = exp (X * b);
          epsilon = 1e-6; % A small constant
          dev = sum(2 * (y .* log((y + epsilon) ./ (p + epsilon)) - (y - p)));
        case "binomial"
          linear_predictor = X * b;
          p = exp(linear_predictor) ./ (1 + exp(linear_predictor)); % Predicted probability of success
          epsilon = 1e-6; % Small constant to prevent log(0)
          if size (y, 2) == 1
            trials = ones(size (y));
          elseif size (y, 2) == 2
            trials = y(:, 2);
          endif
          dev = ...
            -2 * sum(y(:,1) .* log(max(p, epsilon)) + (trials - y(:,1)) .* log(max(1 - p, epsilon)));
        case "normal"
          dev = sum ((y - (X * b)) .^ 2);
      endswitch
      varargout{1} = dev;
  endif
endfunction

%!demo
%! rand("seed", 42); 
%! X = rand(100, 1); 
%! b_true = [0.5; -1.2]; 
%! mu = exp(beta_true(1) + beta_true(2) * X); 
%! y = poissrnd(mu); 
%!
%! ## Fit a GLM model using the Poisson distribution
%! b = glmfit(X, y, "poisson");
%!
%! ## Display the estimated coefficients
%! disp("Estimated coefficients (Intercept, Slope):");
%! disp(b);
%!
%! ## Visualize the fitted model against the original data
%! X_plot = linspace(min(X), max(X), 100)';
%! mu_est = exp([ones(size(X_plot, 1), 1), X_plot] * b); 
%!
%! figure;
%! ## Original data points
%! plot(X, y, 'o', 'MarkerFaceColor', 'w'); hold on; 
%! ## Fitted model curve
%! plot(X_plot, mu_est, 'r-', 'LineWidth', 2); 
%! xlabel('Predictor (X)');
%! ylabel('Response (y)');
%! title('GLM Poisson Regression Fit');
%! legend({'Data', 'Fitted model'}, 'Location', 'NorthEast');
%! hold off;


## Test output
%!test
%! rand ("seed", 1);
%! X = rand (50, 1);
%! b_true = [0.4; 1.5];
%! mu_true = exp (b_true(1) + b_true(2) * X);
%! randp ("seed", 1);
%! y = poissrnd (mu_true);
%! b = glmfit(X, y, "poisson", "link", "log");
%! assert(b(1), b_true(1), 0.5);
%! assert(b(2), b_true(2), 0.5);

## Test input validation
%!error glmfit()
%!error glmfit(rand(5,2))
%!error glmfit(rand(5,2),rand(5,1))
%!error glmfit(rand(5,2), rand(5,1), 'poisson', 'link')
%!error <glmfit: X must be a numeric.> glmfit('abc', rand(6,1), 'poisson')
%!error <glmfit: y must be a numeric.> glmfit(rand(5,2), 'abc', 'poisson')
%!error <glmfit: invalid distribution.> glmfit(rand(5,2), rand(5,1), 2)
%!error <glmfit: invalid parameter name.> glmfit(rand(5,2), rand(5,1), 'poisson', 2, 'log')
%!error <glmfit: invalid link function.> glmfit(rand(5,2), rand(5,1), 'poisson', 'link', 2)
%!error <glmfit: X and y must have same number of observations.> glmfit(rand(5,2), rand(6,1), 'poisson')
%!error <glmfit: unsupported link function.> glmfit(rand(10,2), rand(10,1), 'poisson', 'link', 'inverse')
%!error <glmfit: unsupported parameter name.> glmfit(rand(10,2), rand(10,1), 'poisson', 'notALink', 'log')
%!error <glmfit: distribution must be a recognized value.> glmfit(rand(10,2), rand(10,1), 'xyz')
