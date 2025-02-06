## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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
## @item @var{y} is a @math{n x 1} numeric vector of responses for all supported
## distributions, except for the 'binomial' distribution which can also have
## @var{y} as a @math{n x 2} matrix, where the first column contains the number
## of successes and the second column contains the number of trials.
## @item @var{distribution} specifies the distribution of the response variable
## (e.g., 'poisson').
## @end itemize
##
## @code{@var{b} = glmfit (@dots{}, @var{Name}, @var{Value})}
## specifies additional options using @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{"link"} @tab @tab A character vector specifying a link function
## or a numeric scalar for the 'p' link in the case of the Inverse Gaussian
## distribution. Supported link functions include 'identity', 'log', 'logit',
## 'probit', 'loglog', 'comploglog', 'reciprocal', and 'p'. For custom link
## functions, provide a cell array with three function handles: the link
## function, its derivative, and its inverse.
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
## Supported distributions include 'poisson', 'binomial', 'normal', 'gamma', and
## 'inverse gaussian'. For the inverse Gaussian distribution, the link function
## can be specified as a numeric scalar for the 'p' link.
## @end deftypefn

function [b,varargout] = glmfit (X, y, distribution, varargin)

  ## Check input
  if (nargin < 3)
    error ("glmfit: too few input arguments.");
  elseif (mod (nargin - 3, 2) != 0)
    error ("glmfit: Name-Value arguments must be in pairs.");
  elseif (! isnumeric (X))
    error ("glmfit: X must be a numeric.");
  elseif (! isnumeric (y))
    error ("glmfit: Y must be a numeric.");
  elseif (size (X, 1) != size (y, 1))
    error ("glmfit: X and Y must have the same number of observations.");
  elseif (! ischar (distribution))
    error ("glmfit: DISTRIBUTION must be a character vector.");
  endif

  ## Remove missing values
  xymissing = isnan (y) | any (isnan (X), 2);
  y(xymissing) = [];
  X(xymissing,:) = [];
  [my, ny] = size (y);
  [mx, nx] = size (X);

  ## Check dimensions based on distribution
  if (strcmpi (distribution, 'binomial'))
    if (size (y, 2) > 2)
      error (["glmfit: for a binomial distribution,", ...
              " Y must be an n-by-1 or n-by-2 matrix."]);
    endif
  else
    if (size (y, 2) != 1)
      error (["glmfit: for distributions other than the binomial,", ...
              " Y must be an n-by-1 column vector."]);
    endif
  endif

  ## Add default link functions
  switch (tolower (distribution))
    case "poisson"
      link = "log";
    case "binomial"
      link = "logit";
    case "normal"
      link = "identity";
    case "gamma"
      link = "reciprocal";
    case "inverse gaussian"
      link = -2;
    otherwise
      error ("glmfit: unknown distribution.");
  endswitch

  ## Set default for constant
  constant = "on";

  ## Parse extra parameters
  while (numel (varargin) > 0)
    switch (tolower (varargin {1}))
      case "link"
        linkInput = varargin {2};
        ## Check custom link
        if (iscell (linkInput))
          if (numel (linkInput) != 3)
            error (["glmfit: custom link functions must be", ...
                    " in a three-element cell array."])
          endif
          linkFunc = linkInput{1};
          derLinkFunc = linkInput{2};
          invLinkFunc = linkInput{3};
          ## Check for function_handle
          if (! all (cellfun (@(f) isa (f, 'function_handle'), linkInput)))
            error ("glmfit: custom link functions must be function handles.");
          endif
          ## Test the custom function with a small vector
          try
            testInput = [1; 2; 3; 4; 5];
            testOutput = invLinkFunc (testInput);
            if (! isequal (size (testInput), size (testOutput)))
              error (["glmfit: custom inverse link function must", ...
                      " return output of the same size as input."]);
            endif
          catch
            error (["glmfit: custom inverse link function must", ...
                      " return output of the same size as input."]);
          end_try_catch
          link = "custom";
        ## Check link
        elseif ischar (linkInput) || isstring (linkInput)
          link = tolower (linkInput);
          if (! any (strcmpi (link, {"identity", "log", "logit", "probit", ...
                                     "loglog", "comploglog", "reciprocal", "p"})))
            error ("glmfit: unsupported link function.");
          endif
        elseif isnumeric (linkInput)
          link = "p";
          p_input = linkInput;
        else
          error ("glmfit: invalid value for link function.");
        endif
      case "constant"
        constant = tolower (varargin {2});
        if (! any (strcmpi (constant, {"on", "off"})))
          error ("glmfit: constant should be either 'on' or 'off'.");
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
      ilink = @(x) x;
    case "log"
      ilink = @(x) exp (x);
    case "logit"
      ilink = @(x) exp (x) ./ (1 + exp (x));
    case "probit"
      ilink = @(x) normcdf (x);
    case "loglog"
      ilink = @(x) exp (-exp (x));
    case "comploglog"
      ilink = @(x) 1 - exp (-exp (x));
    case "reciprocal"
      ilink = @(x) 1 ./ x;
    case "p"
      if (isnumeric (p_input))
        ilink = @(x) x .^ p_input;
      else
        error ("glmfit: invalid value for link function.");
      endif
    case "custom"
      ilink = invLinkFunc;
    otherwise
      error ("glmfit: unsupported link function.");
  endswitch

  ## Select negative loglikelihood according to distribution
  switch (tolower (distribution))
    case "poisson"
      nll = @(b) - sum ((y .* X * b) - ilink (X * b) - gammaln (y+1));
    case "binomial"
      eps = 1e-6;
      if (size (y, 2) == 1)
        successes = y;
        trials = ones (size (y));
      else  # it can only have 2 columns then
        successes = y(:, 1);
        trials = y(:, 2);
      endif
      nll = @(b) ...
          - sum (successes .* log (max (min (ilink (X * b), 1 - eps), eps)) ...
          + (trials - successes) ...
          .* log (1 - max (min (ilink (X * b), 1 - eps), eps)));
    case "normal"
      nll = @(b) 0.5 * sum ((y - ilink (X * b)) .^ 2);
    case "gamma"
      nll = @(b) sum ((y ./ ilink (X * b)) + log (ilink (X * b)));
    case "inverse gaussian"
      nll = @(b) sum ((y - ilink (X * b)) .^ 2 ./ (y .* ilink (X * b) .^ 2));
  endswitch

  options = optimset ('MaxFunEvals', 10000, 'MaxIter', 10000);
  b = fminsearch (nll, b, options);

  if (nargout > 1)
    dev = [];
    switch (tolower (distribution))
      case "poisson"
        p = exp (X * b);
        eps = 1e-6;
        dev = sum (2 * (y .* log ((y + eps) ./ (p + eps)) - (y - p)));
      case "binomial"
        eps = 1e-10;
        if (size (y, 2) == 1)
          successes = y;
          trials = ones(size (y));
        elseif (size (y, 2) == 2)
          successes = y(:, 1);
          trials = y(:, 2);
        endif
        p_hat = max (min (ilink (X * b), 1 - eps), eps);
        p = successes ./ trials;
        p = max (min (p, 1 - eps), eps);
        dev = 2 * sum (successes .* log (p ./ p_hat) ...
              + (trials - successes) .* log ((1 - p) ./ (1 - p_hat)));
      case "normal"
        dev = sum ((y - (X * b)) .^ 2);
      case "gamma"
        dev = 2 * sum ((y - ilink (X * b)) ./ ilink (X * b) - log (y ./ ilink (X * b)));
      case "inverse gaussian"
        dev = sum ((y - ilink (X * b)) .^ 2 ./ (y .* ilink (X * b) .^ 2));
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
%!error <glmfit: too few input arguments.> glmfit ()
%!error <glmfit: too few input arguments.> glmfit (1)
%!error <glmfit: too few input arguments.> glmfit (1, 2)
%!error <glmfit: Name-Value arguments must be in pairs.> ...
%! glmfit (rand (6, 1), rand (6, 1), 'poisson', 'link')
%!error <glmfit: X must be a numeric.> ...
%! glmfit ('abc', rand (6, 1), 'poisson')
%!error <glmfit: Y must be a numeric.> ...
%! glmfit (rand (5, 2), 'abc', 'poisson')
%!error <glmfit: X and Y must have the same number of observations.> ...
%! glmfit (rand (5, 2), rand (6, 1), 'poisson')
%!error <glmfit: DISTRIBUTION must be a character vector.> ...
%! glmfit (rand (6, 2), rand (6, 1), 3)
%!error <glmfit: DISTRIBUTION must be a character vector.> ...
%! glmfit (rand (6, 2), rand (6, 1), {'poisson'})
%!error <glmfit: for a binomial distribution, Y must be an n-by-1 or n-by-2 matrix.> ...
%! glmfit (rand (5, 2), rand (5, 3), 'binomial')
%!error <glmfit: for distributions other than the binomial, Y must be an n-by-1 column vector> ...
%! glmfit (rand (5, 2), rand (5, 2), 'normal')
%!error <glmfit: 'gamma' distribution is not supported yet.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'gamma')
%!error <glmfit: 'inverse gaussian' distribution is not supported yet.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'inverse gaussian')
%!error <glmfit: unknown distribution.> ...
%! glmfit (rand (5, 2), rand (5, 1), 'loguniform')
%!error <glmfit: custom link functions must be in a three-element cell array.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {'log'})
%!error <glmfit: custom link functions must be in a three-element cell array.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {'log', 'hijy'})
%!error <glmfit: custom link functions must be function handles.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {'log','dfv','dfgvd'})
%!error <glmfit: custom link functions must be function handles.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {@log, 'derivative', @exp})
%!error <glmfit: custom inverse link function must return output of the same size as input.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', {@exp, @log, @(x) eye(e)})
%!error <glmfit: unsupported link function.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', 'somelinkfunction')
%!error <glmfit: invalid value for link function.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', 2)
%!error <glmfit: constant should be either 'on' or 'off'.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', 'log', 'constant', 0)
%!error <glmfit: constant should be either 'on' or 'off'.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', 'log', 'constant', 'asda')
%!error <glmfit: unknown parameter name.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'param', 'log', 'constant', 'on')
%!error <glmfit: Name-Value arguments must be in pairs.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', 2, 'constant')
%!error <glmfit: constant should be either 'on' or 'off'.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', 2, 'constant', 'invalid')
%!error <glmfit: unknown distribution.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', 2, 'constant', 'invalid')
%!error <glmfit: unsupported link function.> ...
%! glmfit (rand(5,2), rand(5,1), 'poisson', 'link', "2")

%!test
%! rand ("seed", 1);
%! X1 = rand (50, 1);
%! X2 = rand (50, 1) * 0.5;
%! b_true = [0.4; 1.5; -0.7];
%! mu_true = exp (b_true(1) + b_true(2) * X1 + b_true(3) * X2);
%! shape = 2;
%! scale = mu_true ./ shape;
%! y = gamrnd (shape, scale);
%! [b, dev] = glmfit ([X1, X2], y, "gamma", "link", "log");
%! assert (b(1), b_true(1), 0.5);
%! assert (b(2), b_true(2), 0.5);
%! assert (b(3), b_true(3), 0.5);
%! assert (dev < 100, true);

%!test
%! rand ("seed", 1);
%! X1 = rand (50, 1);
%! X2 = rand (50, 1) * 0.5;
%! b_true = [0.4; 1.5; -0.7];
%! mu_true = exp (b_true(1) + b_true(2) * X1 + b_true(3) * X2);
%! lambda = 1;
%! y = invgrnd (mu_true, lambda);
%! [b, dev] = glmfit ([X1, X2], y, "inverse gaussian", "link", "log");
%! assert (b(1), b_true(1), 1.0);
%! assert (b(2), b_true(2), 1.0);
%! assert (b(3), b_true(3), 1.0);
%! assert (dev < 100, true);

%!test
%! rand ("seed", 1);
%! X = rand (50, 1);
%! b_true = [0.4; 1.5];
%! p_input = 2;
%! mu_true = (b_true(1) + b_true(2) * X).^p_input;
%! randp ("seed", 1);
%! y = poissrnd (mu_true);
%! [b, dev] = glmfit (X, y, "poisson", "link", p_input);
%! assert (b(1), b_true(1), 0.7);
%! assert (b(2), b_true(2), 0.7);
%! assert (dev < 100, true);

%!test
%! X = [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 9.0, 10.1]';
%! y = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]';
%! [b, dev] = glmfit (X, y, "gamma", "link", "log");
%! b_matlab = [-0.7631; 0.1113];
%! dev_matlab = 0.0111;
%! assert (b, b_matlab, 0.001);
%! assert (dev, dev_matlab, 0.001);

%!test
%! X = [1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 9.0, 10.1]';
%! y = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]';
%! p_input = 1;
%! [b, dev] = glmfit (X, y, "inverse gaussian", "link", p_input);
%! b_matlab = [0.3813; 0.0950];
%! dev_matlab = 0.0051;
%! assert (b, b_matlab, 0.001);
%! assert (dev, dev_matlab, 0.001);