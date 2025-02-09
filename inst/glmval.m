## Copyright (C) 2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{yhat} =} glmval (@var{b}, @var{X}, @var{link})
## @deftypefnx {statistics} {[@var{yhat}, @var{y_lo}, @var{y_hi}] =} glmval (@var{b}, @var{X}, @var{link}, @var{stats})
## @deftypefnx {statistics} {[@dots{}] =} glmval (@dots{}, @var{Name}, @var{Value})
##
## Predict values for a generalized linear model.
##
## @code{@var{yhat} = glmval (@var{b}, @var{X}, @var{link})} returns the
## predicted values for the generalized linear model with a vector of
## coefficient estimates @var{b}, a matrix of predictors @var{X}, in which each
## column corresponds to a distinct predictor variable, and a link function
## @var{link}, which can be any of the character vectors, numeric scalar, or
## custom-defined link functions used as values for the @qcode{"link"}
## name-value pair argument in the @code{glmfit} function.
##
## @code{[@var{yhat}, @var{y_lo}, @var{y_hi}] = glmval (@var{b}, @var{X},
## @var{link}, @var{stats})} also returns the 95% confidence intervals for the
## predicted values according to the model's statistics contained in the
## @var{stats} structure, which is the output of the @code{glmfit} function.
## By default, the confidence intervals are nonsimultaneous, and apply to the
## fitted curve instead of new observations.
##
## @code{[@dots{}] = glmval (@dots{}, @var{Name}, @var{Value})} specifies
## additional options using @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{"confidence"} @tab @tab A scalar value between 0 and 1
## specifying the confidence level for the confidence bounds.
##
## @item @qcode{"Constant"} @tab @tab A character vector specifying whether to
## include a constant term in the model.  Valid options are @var{"on"} (default)
## and @var{"off"}.
##
## @item @qcode{"simultaneous"} @tab @tab Specifies whether to
## include a constant term in the model. Options are
## @var{"on"} (default) or @var{"off"}.
##
## @item @qcode{"size"} @tab @tab A numeric scalar or a vector with one value
## for each row of @var{X} specifying the size parameter @math{N} for a binomial
## model.
## @end multitable
##
## @seealso{glmfit}
## @end deftypefn

function [yhat, y_lo, y_hi] = glmval (b, X, link, varargin)

  ## Check input arguments
  if (nargin < 3)
    error ("glmval: too few input arguments.");
  elseif (! (isnumeric (b) && isvector (b)) || isempty (b))
    error ("glmval: B must be a numeric vector of coefficient estimates.");
  elseif (! isnumeric (X) || isempty (X))
    error ("glmval: X must be a numeric matrix.");
  endif

  ## Get inverse link (input validation is performed in private function)
  [~, ~, ilink, errmsg] = getlinkfunctions (link);
  if (! isempty (errmsg))
    error ("glmval: %s", errmsg);
  endif

  ## Check if fourth input argument is a STATS structure
  stats = [];
  if (nargin > 3)
    if (isstruct (varargin{1}))
      stats = varargin{1};
      rf = {"s", "se", "coeffcorr", "estdisp", "dfe"};
      if (! all (ismember (rf, fieldnames (stats))))
        error ("glmval: invalid 'stats' structure.");
      endif
      varargin(1) = [];
    endif
  endif

  ## Set defaults
  confidence = 0.95;
  constant = true;
  offset = zeros (size (X, 1), 1);
  simultaneous = false;
  N = 1;

  ## Parse extra parameters
  if (mod (numel (varargin), 2) != 0)
    error ("glmval: Name-Value arguments must be in pairs.");
  endif
  while (numel (varargin) > 0)
    switch (tolower (varargin {1}))

      case "confidence"
        confidence = varargin {2};
        if (! (isscalar (confidence) && isnumeric (confidence)
                                     && confidence > 0 && confidence < 1))
          error ("glmval: 'Confidence' must be a scalar between 0 and 1.");
        endif

      case "constant"
        constant = tolower (varargin {2});
        if (strcmpi (constant, "on"))
          constant = true;
        elseif (strcmpi (constant, "off"))
          constant = false;
        else
          error ("glmval: 'Constant' should be either 'on' or 'off'.");
        endif

      case "offset"
        offset = varargin {2};
        if (! (isnumeric (offset) && isequal (numel (offset), size (X, 1))))
          error (["glmval: 'Offset' must be a numeric vector", ...
                  " of the same length as the rows in X."]);
        endif
        offset = offset(:);

      case "simultaneous"
        simultaneous = varargin {2};
        if (! (islogical (simultaneous) && isscalar (simultaneous)))
          error ("glmval: 'simultaneous' must be a boolean scalar.");
        endif

      case "size"
        N = varargin {2};
        if (! isnumeric (N) ||
            ! (isscalar (N) || isvector (N) && isequal (numel (N), size (X, 1))))
          error (["glmval: 'size' must be a scalar or a vector with", ...
                  " one value for each row of X."]);
        endif
        N = N(:);

      otherwise
        error ("glmval: unknown parameter name.");
    endswitch
    varargin (1:2) = [];
  endwhile

  ## Adjust X based on constant
  if (constant)
    X = [ones(size (X, 1), 1), X];
  endif

  ## Predict yhat
  eta = X * b + offset;
  yhat = N .* ilink (eta);

  ## Compute lower and upper bounds
  if (nargout > 1)
    if (isempty (stats))
      error (["glmval: cannot compute confidence", ...
              " intervals without STATS structure."]);
    endif
    if (isnan (stats.s))
      y_lo = NaN (size (yhat));
      y_hi = NaN (size (yhat));
    else
      V = (stats.se * stats.se') .* stats.coeffcorr;
      XVX = sum ((X * V) .* X, 2)';
      if (simultaneous)
        dof = length (b);
        if (stats.estdisp)
          Xcv = sqrt (dof * finv (confidence, dof, stats.dfe));
        else
          Xcv = sqrt (chi2inv (confidence, dof));
        endif
      else
        Xcv = tinv ((1 + confidence) / 2, stats.dfe);
      endif
      interval = Xcv * sqrt (XVX);
      int_hilo = [N.*ilink(eta-interval) N.*ilink(eta+interval)];
      y_lo = yhat - min (int_hilo, [], 2);
      y_hi = max (int_hilo, [], 2) - yhat;
    endif
  endif

endfunction

%!demo
%! x = [210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430]';
%! n = [48, 42, 31, 34, 31, 21, 23, 23, 21, 16, 17, 21]';
%! y = [1, 2, 0, 3, 8, 8, 14, 17, 19, 15, 17, 21]';
%! b = glmfit (x, [y n], "binomial", "Link", "probit");
%! yfit = glmval (b, x, "probit", "Size", n);
%! plot (x, y./n, 'o', x, yfit ./ n, '-')

## Test input validation
%!error <glmval: too few input arguments.> glmval ()
%!error <glmval: too few input arguments.> glmval (1)
%!error <glmval: too few input arguments.> glmval (1, 2)
%!error <glmval: B must be a numeric vector of coefficient estimates.> ...
%! glmval ("asd", [1; 1; 1], 'probit')
%!error <glmval: B must be a numeric vector of coefficient estimates.> ...
%! glmval ([], [1; 1; 1], 'probit')
%!error <glmval: X must be a numeric matrix.> ...
%! glmval ([0.1; 0.3; 0.4], [], 'probit')
%!error <glmval: X must be a numeric matrix.> ...
%! glmval ([0.1; 0.3; 0.4], "asd", 'probit')
%!error <glmval: structure with custom link functions must be a scalar.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", {1, 2}))
%!error <glmval: structure with custom link functions requires the fields 'Link', 'Derivative', and 'Inverse'.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", "norminv"))
%!error <glmval: bad 'Link' function in custom link function structure.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", "some", "Derivative", @(x)x, "Inverse", "normcdf"))
%!error <glmval: bad 'Link' function in custom link function structure.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", 1, "Derivative", @(x)x, "Inverse", "normcdf"))
%!error <glmval: custom 'Link' function must return an output of the same size as input.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", @(x) [x, x], "Derivative", @(x)x, "Inverse", "normcdf"))
%!error <glmval: invalid custom 'Link' function.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", "what", "Derivative", @(x)x, "Inverse", "normcdf"))
%!error <glmval: bad 'Derivative' function in custom link function structure.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", @(x)x, "Derivative", "some", "Inverse", "normcdf"))
%!error <glmval: bad 'Derivative' function in custom link function structure.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", @(x)x, "Derivative", 1, "Inverse", "normcdf"))
%!error <glmval: custom 'Derivative' function must return an output of the same size as input.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", @(x)x, "Derivative", @(x) [x, x], "Inverse", "normcdf"))
%!error <glmval: invalid custom 'Derivative' function.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", @(x)x, "Derivative", "what", "Inverse", "normcdf"))
%!error <glmval: bad 'Inverse' function in custom link function structure.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", @(x)x, "Derivative", "normcdf", "Inverse", "some"))
%!error <glmval: bad 'Inverse' function in custom link function structure.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", @(x)x, "Derivative", "normcdf", "Inverse", 1))
%!error <glmval: custom 'Inverse' function must return an output of the same size as input.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", @(x)x, "Derivative", "normcdf", "Inverse", @(x) [x, x]))
%!error <glmval: invalid custom 'Inverse' function.> ...
%! glmval (rand (3,1), rand (5,2), struct ("Link", @(x)x, "Derivative", "normcdf", "Inverse", "what"))
%!error <glmval: cell array with custom link functions must have three elements.> ...
%! glmval (rand (3,1), rand (5,2), {'log'})
%!error <glmval: cell array with custom link functions must have three elements.> ...
%! glmval (rand (3,1), rand (5,2), {'log', 'hijy'})
%!error <glmval: cell array with custom link functions must have three elements.> ...
%! glmval (rand (3,1), rand (5,2), {1, 2, 3, 4})
%!error <glmval: bad 'Link' function in custom link function cell array.> ...
%! glmval (rand (3,1), rand (5,2), {"log", "dfv", "dfgvd"})
%!error <glmval: custom 'Link' function must return an output of the same size as input.> ...
%! glmval (rand (3,1), rand (5,2), {@(x) [x, x], "dfv", "dfgvd"})
%!error <glmval: invalid custom 'Link' function.> ...
%! glmval (rand (3,1), rand (5,2), {@(x) what (x), "dfv", "dfgvd"})
%!error <glmval: bad 'Derivative' function in custom link function cell array.> ...
%! glmval (rand (3,1), rand (5,2), {@(x) x, "dfv", "dfgvd"})
%!error <glmval: custom 'Derivative' function must return an output of the same size as input.> ...
%! glmval (rand (3,1), rand (5,2), {@(x) x, @(x) [x, x], "dfgvd"})
%!error <glmval: invalid custom 'Derivative' function.> ...
%! glmval (rand (3,1), rand (5,2), {@(x) x, @(x) what (x), "dfgvd"})
%!error <glmval: bad 'Inverse' function in custom link function cell array.> ...
%! glmval (rand (3,1), rand (5,2), {@(x) x, @(x) x, "dfgvd"})
%!error <glmval: custom 'Inverse' function must return an output of the same size as input.> ...
%! glmval (rand (3,1), rand (5,2), {@(x) x, @(x) x, @(x) [x, x]})
%!error <glmval: invalid custom 'Inverse' function.> ...
%! glmval (rand (3,1), rand (5,2), {@(x) x, @(x) x, @(x) what (x)})
%!error <glmval: numeric input for custom link function must be a finite real scalar value.> ...
%! glmval (rand (3,1), rand (5,2), NaN)
%!error <glmval: numeric input for custom link function must be a finite real scalar value.> ...
%! glmval (rand (3,1), rand (5,2), [1, 2])
%!error <glmval: numeric input for custom link function must be a finite real scalar value.> ...
%! glmval (rand (3,1), rand (5,2), [1i])
%!error <glmval: canonical link function name must be a character vector.> ...
%! glmval (rand (3,1), rand (5,2), ["log"; "log1"])
%!error <glmval: canonical link function 'somelinkfunction' is not supported.> ...
%! glmval (rand (3,1), rand (5,2), 'somelinkfunction')
%!error <glmval: invalid value for custom link function.> ...
%! glmval (rand (3,1), rand (5,2), true)
%!error <glmval: invalid 'stats' structure.> ...
%! glmval (rand (3,1), rand (5,2), 'probit', struct ("s", 1))
%!error <glmval: Name-Value arguments must be in pairs.> ...
%! glmval (rand (3,1), rand (5,2), 'probit', 'confidence')
%!error <glmval: 'Confidence' must be a scalar between 0 and 1.> ...
%! glmval (rand (3,1), rand (5,2), 'probit', 'confidence', 0)
%!error <glmval: 'Confidence' must be a scalar between 0 and 1.> ...
%! glmval (rand (3,1), rand (5,2), 'probit', 'confidence', 1.2)
%!error <glmval: 'Confidence' must be a scalar between 0 and 1.> ...
%! glmval (rand (3,1), rand (5,2), 'probit', 'confidence', [0.9, 0.95])
%!error <glmval: 'Constant' should be either 'on' or 'off'.> ...
%! glmval (rand (3, 1), rand (5, 2), 'probit', 'constant', 1)
%!error <glmval: 'Constant' should be either 'on' or 'off'.> ...
%! glmval (rand (3, 1), rand (5, 2), 'probit', 'constant', 'o')
%!error <glmval: 'Constant' should be either 'on' or 'off'.> ...
%! glmval (rand (3, 1), rand (5, 2), 'probit', 'constant', true)
%!error <glmval: 'Offset' must be a numeric vector of the same length as the rows in X.> ...
%! glmval (rand (3, 1), rand (5, 2), 'probit', 'offset', [1; 2; 3; 4])
%!error <glmval: 'Offset' must be a numeric vector of the same length as the rows in X.> ...
%! glmval (rand (3, 1), rand (5, 2), 'probit', 'offset', 'asdfg')
%!error <glmval: 'simultaneous' must be a boolean scalar.> ...
%! glmval (rand (3, 1), rand (5, 2), 'probit', 'simultaneous', 'asdfg')
%!error <glmval: 'simultaneous' must be a boolean scalar.> ...
%! glmval (rand (3, 1), rand (5, 2), 'probit', 'simultaneous', [true, false])
%!error <glmval: 'size' must be a scalar or a vector with one value for each row of X.> ...
%! glmval (rand (3, 1), rand (5, 2), 'probit', 'size', "asd")
%!error <glmval: 'size' must be a scalar or a vector with one value for each row of X.> ...
%! glmval (rand (3, 1), rand (5, 2), 'probit', 'size', [2, 3, 4])
%!error <glmval: 'size' must be a scalar or a vector with one value for each row of X.> ...
%! glmval (rand (3, 1), rand (5, 2), 'probit', 'size', [2; 3; 4])
%!error <glmval: 'size' must be a scalar or a vector with one value for each row of X.> ...
%! glmval (rand (3, 1), rand (5, 2), 'probit', 'size', ones (3))
%!error <glmval: unknown parameter name.> ...
%! glmval (rand (3, 1), rand (5, 2), 'probit', 'someparam', 4)
%!error <glmval: cannot compute confidence intervals without STATS structure.> ...
%! [y,lo,hi] = glmval (rand (3, 1), rand (5, 2), 'probit')
