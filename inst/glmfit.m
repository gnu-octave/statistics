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
## @deftypefn {statistics} {@var{b} =} glmfit (@var{X}, @var{y}, @var{distribution})
## @deftypefnx {statistics} {@var{b} =} glmfit (@var{X}, @var{y}, @var{distribution},@var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{b}, @var{dev}] =} glmfit (@dots{})
##
## Perform generalized linear model fitting.
##
## @code{@var{b} = glmfit (@var{X}, @var{y}, @var{distribution})} returns a 
## coefficient estimates vector, @var{b} for a
## generalized linear regression model of responses in @var{y} and 
## predictors in @var{X}, using the @var{distribution}.
##
## @itemize
## @item @var{X} is an @math{nxp} numeric matrix of predictor variables with
## @math{n} observations and @math{p} predictors.
## @item @var{y} is a @math{n x 1} vector of responses for all the distributions. 
## 'binomial' distribution can have @var{y} as a @math{n x 2} matrix 
## where the first column contains the number of successes and the
## second column contains the number of trials.
## @item @var{distribution} specifies the distribution of the response variable
## (e.g., 'poisson').
## @end itemize
##
## @code{@var{b} = glmfit (@dots{}, @var{Name}, @var{Value})} 
## specifies additional options using @qcode{Name-Value} pair arguments.
## 
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @qcode{"link"} @tab @tab A character vector specifying a link
## function.  
##
## @item @qcode{"constant"} @tab @tab Specifies whether to 
## include a constant term in the model. Options are 
## @var{"on"} (default) or @var{"off"}.
## @end multitable
##
## @code{[@var{b}, @var{dev}] = glmfit (@dots{})} 
## returns the estimated coefficient vector, @var{b}, as well as 
## the deviance, @var{dev}, of the fit.
##
## Supported distributions include 'poisson', 'binomial', and 'normal'. 
## Supported link functions include 'identity', 'log', 'logit', 'probit', 
## 'loglog', 'comploglog', 'reciprocal' and a custom link.
## Custom link function provided as a structure with three fields:
## Link Function, Derivative Function, Inverse Function. 
## @end deftypefn

function [b,varargout] = glmfit (X, y, distribution, varargin)
  ## Check input
  if nargin < 3
    error ("glmfit: at least three input arguments (X, y, distribution) are required.");
  elseif mod (nargin - 3, 2) != 0
    error ("glmfit: Name-Value arguments must be in pairs.");
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
  if strcmpi (distribution, 'binomial')
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
          if ! all (cellfun (@(f) isa (f, 'function_handle'), linkInput))
            error ("glmfit: custom link functions must be function handles.");
          endif
          ## Test the custom function with a small vector
          try
            testInput = [1; 2; 3; 4; 5];
            testOutput = invLinkFunc (testInput);
            if ! isequal (size (testInput), size (testOutput))
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
  if (strcmpi (constant, 'on'))
    X = [ones(mx, 1), X];
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
  b = fminsearch (negLogLikelihood, b, options);
  
  if (nargout > 1)
      dev = [];
      switch (tolower (distribution))
        case "poisson"
          p = exp (X * b);
          eps = 1e-6; 
          dev = sum (2 * (y .* log ((y + eps) ./ (p + eps)) - (y - p)));
        case "binomial"
          eps = 1e-10; 
          if size (y, 2) == 1
            successes = y;
            trials = ones(size (y));
          elseif size (y, 2) == 2
            successes = y(:, 1);
            trials = y(:, 2);
          endif
          p_hat = max (min (inverse_link_func (X * b), 1 - eps), eps);
          p = successes ./ trials;
          p = max (min (p, 1 - eps), eps);
          dev = 2 * sum (successes .* log (p ./ p_hat) + (trials - successes) ...
            .* log ((1 - p) ./ (1 - p_hat)));
        case "normal"
          dev = sum ((y - (X * b)) .^ 2);
      endswitch
      varargout{1} = dev;
  endif
endfunction

%!demo
%! rand ("seed", 1); 
%! X = rand (100, 1); 
%! b_true = [0.5; -1.2]; 
%! mu = exp (b_true(1) + b_true(2) * X); 
%! randp ("seed", 1);
%! y = poissrnd (mu); 
%! ## Fit a GLM model using the poisson distribution
%! [b,dev] = glmfit (X, y, 'poisson');

%!demo
%! x = [2100 2300 2500 2700 2900 3100 3300 3500 3700 3900 4100 4300]';
%! n = [48 42 31 34 31 21 23 23 21 16 17 21]';
%! y = [1 2 0 3 8 8 14 17 19 15 17 21]';
%! [b,dev] = glmfit (x,[y n],'binomial','Link','probit');

## Test output
%!test
%! rand ("seed", 1);
%! X = rand (50, 1);
%! b_true = [0.4; 1.5];
%! mu_true = exp (b_true(1) + b_true(2) * X);
%! randp ("seed", 1);
%! y = poissrnd (mu_true);
%! b = glmfit (X, y, "poisson", "link", "log");
%! assert (b(1), b_true(1), 0.5);
%! assert (b(2), b_true(2), 0.5);
%!test
%! rand ("seed", 1);
%! X1 = rand (50, 1);
%! X2 = rand (50, 1) * 0.5; 
%! b_true = [0.4; 1.5; -0.7]; 
%! mu_true = exp (b_true(1) + b_true(2) * X1 + b_true(3) * X2);
%! randp ("seed", 1);
%! y = poissrnd(mu_true);
%! [b, dev] = glmfit ([X1, X2], y, "poisson", "link", "log");
%! assert (b(1), b_true(1), 1);
%! assert (b(2), b_true(2), 1);
%! assert (b(3), b_true(3), 1);
%! assert (dev < 60, true); 

## Test input validation
%!error <glmfit: too few arguments.> ...
%! glmfit ()
%!error <glmfit: too few arguments.> ...
%! glmfit (rand(5,2))
%!error <glmfit: too few arguments.> ...
%! glmfit (rand(5,2),rand(5,1))
%!error <glmfit: too few arguments.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link')
%!error <glmfit: X must be a numeric.> ...
%! glmfit ('abc', rand(6,1), 'poisson')
%!error <glmfit: y must be a numeric.> ...
%! glmfit (rand(5,2), 'abc', 'poisson')
%!error <glmfit: invalid distribution.> ...
%! glmfit (rand(5,2), rand(5,1), 2)
%!error <glmfit: invalid parameter name.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 2, 'log')
%!error <glmfit: invalid link function.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', 2)
%!error <glmfit: X and y must have same number of observations.> ...
%! glmfit (rand(5,2), rand(6,1), 'poisson')
%!error <glmfit: unsupported link function.> ...
%! glmfit (rand(10,2), rand(10,1), 'poisson', 'link', 'inverse')
%!error <glmfit: unsupported parameter name.> ...
%! glmfit (rand(10,2), rand(10,1), 'poisson', 'notALink', 'log')
%!error <glmfit: distribution must be a recognized value.> ...
%! glmfit (rand(10,2), rand(10,1), 'xyz')
%!error <glmfit: y must be a vector or n-by-2 matrix for binomial distribution.> ...
%! glmfit (rand(5,2), rand(5,3), 'binomial')
%!error <glmfit: constant value should be either 'on' or 'off'.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'constant', 'yes')
%!error <glmfit: y cannot contain negative values for 'binomial' distribution.> ...
%! glmfit (rand(5,2), [-1, 2; 2, 3; 3, 4; 4, 5; 5, 6], 'binomial')
%!error <glmfit: number of successes cannot exceed number of trials.> ...
%! glmfit (rand(5,2), [randi(10,5,1), randi(9,5,1)], 'binomial')
%!error <glmfit: y cannot contain negative values for 'poisson' distribution.> ...
%! glmfit (rand(5,2), [-1; 2; 3; 4; 5], 'poisson')
%!error <glmfit: custom link function must include three function handles.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {'log'})
%!error <glmfit: custom link function must include three function handles.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {'log', 'hijy'})
%!error <glmfit: custom link function must include three function handles.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {'log','dfv','dfgvd'})
%!error <glmfit: for custom link, all inputs must be function handles.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {@log, 'derivative', @exp})
%!error <glmfit: for custom link, all inputs must be function handles.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {'derivative', @log, @exp})
%!error <glmfit: for custom link, all inputs must be function handles.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {@log, @exp, 'derivative'})
