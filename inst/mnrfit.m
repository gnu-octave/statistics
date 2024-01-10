## Copyright (C) 2024 Andrew C Penn
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
## @deftypefn  {statistics} {@var{B} =} mnrfit (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{B} =} mnrfit (@var{X}, @var{Y}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{B}, @var{dev}] =} anovan (@dots{})
## @deftypefnx {statistics} {[@var{B}, @var{dev}, @var{stats}] =} anovan (@dots{})
##
## Perform logistic regression for binomial responses or multiple ordinal
## responses.
##
## @seealso{logistic_regression}
## @end deftypefn


function [B, DEV, STATS] = mnrfit (X, Y, varargin)

  ## Check supplied parameters
  if ((numel (varargin) / 2) != fix (numel (varargin) / 2))
    error ("mnrfit: wrong number of arguments.")
  endif
  MODELTYPE = "nominal";
  DISPLAY = "off";
  for idx = 3:2:nargin
    name = varargin{idx-2};
    value = varargin{idx-1};
    switch (lower (name))
      case "model"
        MODELTYPE = value;
      case "display"
        DISPLAY = value;
      otherwise
        warning (sprintf ("mnrfit: parameter %s will be ignored", name));
      endswitch
    endfor

  ## Evaluate display input argument
  switch (lower (DISPLAY))
    case "on"
      dispopt = true;
    case "off"
      dispopt = false;
  end

  ## Evaluate X
  if (! (isa (X, 'single') || isa (X, 'double')))
    error ("mnrfit: Predictors must be numeric.")
  endif
  if (isscalar (X) || (ndims (X) > 2))
    error ("mnrfit: Predictors must be a vector or 2D-matrix.")
  endif
  [N, P] = size (X);
  if (N == 1)
    % If X is a row vector, make it a column vector
    X = X.';
    N = P;
    P = 1;
  endif

  ## Evaluate Y
  if (! (isnumeric (Y) || islogical (Y) || ischar (Y) || iscellstr (Y)))
    error (sprintf (strcat ("mnrfit: Response labels must be a character", ...
                    " array, a cell array of strings, \nor a vector or", ...
                    " matrix of doubles, singles or logical values.")));
  endif
  if (isscalar (Y) || (ndims (Y) > 2))
    error ("mnrfit: Response must be a vector or 2D-matrix.")
  endif
  if (size (Y, 1) ~= N)
    error ("mnrfit: Y must have the same number of rows as X")
  endif
  K = size (Y, 2);
  if (K > 1)
    ## if K > 1, we are expecting Y to be a logical matrix, although we will
    ## tolerate the equivalent in singles or doubles
    if (! (isnumeric (Y) || islogical (Y)))
      error (sprintf (strcat ("mnrfit: If K > 2, Y must be a vector or", ...
                    " matrix of doubles, singles or \nlogical values")))
    endif
    if (! (all (double(Y(:)) == Y(:))))
      error ("mnrfit: If Y is a matrix, it should only contain 0 and 1")
    endif
    ## Convert Y to a vector of positive integer categories
    Y = sum (bsxfun (@times, (1:K), Y), 2);
  endif
  ridx = any ([any(isnan(X),2), any(isnan(Y),2)], 2);
  Y(ridx, :) = [];           ## remove NaN
  X(ridx, :) = [];           ## remove NaN
  [UY, ~, YN] = unique (Y);  ## find unique categories in the response
  n = numel (UY);            ## number of unique response categories
  if (isnumeric (Y))
    if (! all (Y == YN))
      error ("mnrfit: Y must contain positive integer category numbers")
    endif
  endif

  ## Evaluate model type input argument
  switch (lower (MODELTYPE))
    case "nominal"
      if (n > 2)
        error ("mnrfit: Fitting more than 2 nominal responses are not suppored")
      else
        ## Y has two responses. Ordinal logistic regression can be used to fit
        ## models with binary nominal responses
      endif
    case "ordinal"
      ## Do nothing, ordinal responses are fully supported
    case "hierarchical"
      error ("mnrfit: Fitting hierarchical responses are not suppored")
    otherwise
      error ("mnrfit: Model type not recognised")
  endswitch

  ## Perform fit and reformat output
  [INTERCEPT, SLOPE, DEV, ~, ~, ~, S] = logistic_regression (YN - 1, X, dispopt);
  B = cat (1, INTERCEPT, SLOPE);
  STATS = struct ("beta", B, ...
                  "dfe", [], ...      ## Not used
                  "s", [], ...        ## Not used
                  "sfit", [], ...     ## Not used
                  "estdisp", [], ...  ## Not used
                  "covb", S.cov, ...
                  "coeffcorr", S.coeffcorr, ...
                  "se", S.se, ...
                  "t", [], ...        ## Not used
                  "p", [], ...        ## Not used
                  "resid", [], ...    ## Not used
                  "residp", [], ...   ## Not used
                  "residd", []);      ## Not used

endfunction

