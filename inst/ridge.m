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

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{b} =} ridge (@var{y}, @var{X}, @var{k})
## @deftypefnx {statistics} {@var{b} =} ridge (@var{y}, @var{X}, @var{k}, @var{scaled})
##
## Ridge regression.
##
## @code{@var{b} = ridge (@var{y}, @var{X}, @var{k})} returns the vector of
## coefficient estimates by applying ridge regression from the predictor matrix
## @var{X} to the response vector @var{y}.  Each value of @var{b} is the
## coefficient for the respective ridge parameter given @var{k}.  By default,
## @var{b} is calculated after centering and scaling the predictors to have a
## zero mean and standard deviation 1.
##
## @code{@var{b} = ridge (@var{y}, @var{X}, @var{k}, @var{scaled})} performs the
## regression with the specified scaling of the coefficient estimates @var{b}.
## When @qcode{@var{scaled} = 0}, the function restores the coefficients to the
## scale of the original data thus is more useful for making predictions.  When
## @qcode{@var{scaled} = 1}, the coefficient estimates correspond to the scaled
## centered data.
##
## @itemize
## @item
## @code{y} must be an @math{Nx1} numeric vector with the response data.
## @item
## @code{X} must be an @math{Nxp} numeric matrix with the predictor data.
## @item
## @code{k} must be a numeric vectir with the ridge parameters.
## @item
## @code{scaled} must be a numeric scalar indicating whether the coefficient
## estimates in @var{b} are restored to the scale of the original data.  By
## default, @qcode{@var{scaled} = 1}.
## @end itemize
##
## Further information about Ridge regression can be found at
## @url{https://en.wikipedia.org/wiki/Ridge_regression}
##
## @seealso{lasso, stepwisefit, regress}
## @end deftypefn

function b = ridge (y, X, k, scaled)

  ## Check input arguments
  if (nargin < 3)
    error ("ridge: function called with too few input arguments.");
  endif
  if (! isvector (y) || columns (y) != 1 || isempty (y))
    error ("ridge: Y must be a numeric column vector.");
  endif
  if (! ismatrix (X) || isempty (X))
    error ("ridge: X must be a numeric matrix.");
  endif
  if (rows (y) != rows (X))
    error ("ridge: Y and X must contain the same number of rows.");
  endif

  ## Parse 4th input argument
  if (nargin < 4 || isempty (scaled))
    unscale = false;
  elseif (scaled == 1)
    unscale = false;
  elseif (scaled == 0)
    unscale = true;
  else
    error ("ridge: wrong value for SCALED argument.");
  endif

  ## Force y to a column vector
  y = y(:);

  ## Rremove any missing values
  notnans = ! logical (sum (isnan ([y, X]), 2));
  y = y(notnans);
  X = X(notnans,:);

  ## Scale and center X to zero mean and StD = 1
  m = mean (X);
  stdx = std (X, 0, 1);

  z = (X - m) ./ stdx;

  ## Add pseudo observations
  Z_pseudo = [z; (sqrt(k(1)) .* eye (columns(X)))];
  Y_pseudo = [y; zeros(columns(X), 1)];

  ## Compute coefficients
  b = Z_pseudo \ Y_pseudo;

  nk = numel (k);

  ## Compute the coefficient estimates for additional ridge parameters.
  if (nk >= 2)

    ## Adding a multiple of the identity matrix to the last p rows.

    ## b is set to 0 for the current ridge parameter value
    b(end,nk) = 0;

    for i=2:nk
      Z_pseudo(end-columns(X)+1:end, :)  =  sqrt (k(i)) .* eye (columns (X));
      b(:,i) = Z_pseudo \ Y_pseudo;
    endfor
  endif

  ## Changing back to the scale
  if (unscale)
    b = b ./ repmat (stdx', 1, nk);
    b = [mean(y)-m.*b; b];
  endif

endfunction

%!demo
%! ## Add some description here!!
%!
%! load acetylene
%!
%! X = [x1, x2, x3];
%!
%! x0 = ones (16, 1);
%! x1x2 = x1 .* x2;
%! x1x3 = x1 .* x3;
%! x2x3 = x2 .* x3;
%!
%! D = [x1, x2, x3, x1x2, x1x3, x2x3];
%!
%! k = 0:1e-5:5e-3;
%!
%! b = ridge (y, D, k);
%!
%! figure
%! plot (k, b, "LineWidth", 2)
%! ylim ([-100, 100])
%! grid on
%! xlabel ("Ridge Parameter")
%! ylabel ("Standardized Coefficient")
%! title ("Ridge Trace")
%! legend ("x1", "x2", "x3", "x1x2", "x1x3", "x2x3")

## Test input


## Test input validation
%!error<ridge: function called with too few input arguments.> ridge (1)
%!error<ridge: function called with too few input arguments.> ridge (1, 2)
%!error<ridge: Y must be a numeric column vector.> ridge (ones (3), ones (3), 2)
%!error<ridge: Y must be a numeric column vector.> ridge ([1, 2], ones (2), 2)
%!error<ridge: Y must be a numeric column vector.> ridge ([], ones (3), 2)
%!error<ridge: X must be a numeric matrix.> ridge (ones (5,1), [], 2)
%!error<ridge: Y and X must contain the same number of rows.> ...
%! ridge ([1; 2; 3; 4; 5], ones (3), 3)
%!error<ridge: wrong value for SCALED argument.> ...
%! ridge ([1; 2; 3], ones (3), 3, 2)
%!error<ridge: wrong value for SCALED argument.> ...
%! ridge ([1; 2; 3], ones (3), 3, "some")
