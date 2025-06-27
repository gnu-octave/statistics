## Copyright (C) 2024 Andrew C Penn
## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefnx {statistics} {[@var{B}, @var{dev}] =} mnrfit (@dots{})
## @deftypefnx {statistics} {[@var{B}, @var{dev}, @var{stats}] =} mnrfit (@dots{})
##
## Perform logistic regression for binomial responses or multiple ordinal
## responses.
##
## Note: This function is currently a wrapper for the @code{logistic_regression}
## function. It can only be used for fitting an ordinal logistic model and a
## nominal model with 2 categories (which is an ordinal case).  Hierarchical
## models as well as nominal model with more than two classes are not currently
## supported.  This function is a work in progress.
##
## @code{@var{B} = mnrfit (@var{X}, @var{Y})}  returns a matrix, @var{B}, of
## coefficient estimates for a multinomial logistic regression of the nominal
## responses in @var{Y} on the predictors in @var{X}.  @var{X} is an @math{NxP}
## numeric matrix the observations on predictor variables, where @math{N}
## corresponds to the number of observations and @math{P} corresponds to
## predictor variables.  @var{Y} contains the response category labels and it
## either be an @math{NxP} categorical or numerical matrix (containing only 1s
## and 0s) or an @math{Nx1} numeric vector with positive integer values, a cell
## array of character vectors and a logical vector.  @var{Y} can also be defined
## as a character matrix with each row corresponding to an observation of
## @var{X}.
##
## @code{@var{B} = mnrfit (@var{X}, @var{Y}, @var{name}, @var{value})} returns a
## matrix, @var{B}, of coefficient estimates for a multinomial model fit with
## additional parameterss specified @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{"model"} @tab @tab Specifies the type of model to fit.
## Currently, only @qcode{"ordinal"} is fully supported.  @qcode{"nominal"} is
## only supported for 2 classes in @var{Y}.
##
## @item @qcode{"display"} @tab @tab A flag to enable/disable displaying
## information about the fitted model.  Default is @qcode{"off"}.
## @end multitable
##
## @code{[@var{B}, @var{dev}, @var{stats}] = mnrfit (@dots{}}  also returns the
## deviance of the fit, @var{dev}, and the structure @var{stats} for any of the
## previous input arguments. @var{stats} currently only returns values for the
## fields @qcode{"beta"}, same as @var{B}, @qcode{"coeffcorr"}, the estimated
## correlation matrix for @var{B}, @qcode{"covd"}, the estimated covariance
## matrix for @var{B}, and @qcode{"se"}, the standard errors of the coefficient
## estimates @var{B}.
##
## @seealso{logistic_regression}
## @end deftypefn

function [B, DEV, STATS] = mnrfit (X, Y, varargin)

  ## Check input arguments X and Y
  if (nargin < 2)
    error ("mnrfit: too few input arguments.");
  endif
  if (! isnumeric (X))
    error ("mnrfit: Predictors must be numeric.")
  endif
  if (isscalar (X) || (ndims (X) > 2))
    error ("mnrfit: Predictors must be a vector or a 2D matrix.")
  endif
  if (isscalar (Y) || (ndims (Y) > 2))
    error ("mnrfit: Response must be a vector or a 2D matrix.")
  endif
  [N, P] = size (X);
  [n, K] = size (Y);
  if (N == 1)   ## if X is a row vector, make it a column vector
    X = X(:);
    N = P;
    P = 1;
  endif
  if (n != N)
    error ("mnrfit: Y must have the same number of rows as X.")
  endif
  if (! (isnumeric (Y) || islogical (Y) || ischar (Y) || iscellstr (Y)))
    error (strcat (["mnrfit: Response labels must be a character array,"], ...
                   [" a cell vector of strings, \nor a vector or"], ...
                   [" matrix of doubles, singles or logical values."]));
  endif

  ## Check supplied parameters
  if (mod (numel (varargin), 2) != 0)
    error ("mnrfit: optional arguments must be in pairs.")
  endif
  MODELTYPE = "nominal";
  DISPLAY = "off";
  while (numel (varargin) > 0)
    name = varargin{1};
    value = varargin{2};
    switch (lower (name))
      case "model"
        MODELTYPE = value;
      case "display"
        DISPLAY = value;
      otherwise
        warning (sprintf ("mnrfit: parameter %s will be ignored", name));
    endswitch
    varargin (1:2) = [];
  endwhile

  ## Evaluate display input argument
  switch (lower (DISPLAY))
    case "on"
      dispopt = true;
    case "off"
      dispopt = false;
  endswitch

  ## Categorize Y if it is a cellstring array
  if (iscellstr (Y))
    if (K > 1)
      error ("mnrfit: Y must be a column vector when given as cellstr.");
    endif
    ## Get groups in Y
    [YN, ~, UY] = grp2idx (Y);  # this will also catch "" as missing values
    ## Remove missing values from X and Y
    RowsUsed  = ! logical (sum (isnan ([X, YN]), 2));
    Y         = Y (RowsUsed);
    X         = X (RowsUsed, :);
    ## Renew groups in Y
    [YN, ~, UY] = grp2idx (Y);  # in case a category is removed due to NaNs in X
    n = numel (UY);
  endif

  ## Categorize Y if it is a character array
  if (ischar (Y))
    ## Get groups in Y
    [YN, ~, UY] = grp2idx (Y);  # this will also catch "" as missing values
    ## Remove missing values from X and Y
    RowsUsed  = ! logical (sum (isnan ([X, YN]), 2));
    Y         = Y (RowsUsed);
    X         = X (RowsUsed, :);
    ## Renew groups in Y
    [YN, ~, UY] = grp2idx (Y);  # in case a category is removed due to NaNs in X
    n = numel (UY);
  endif

  if (K > 1)
    ## So far, if K > 1, Y must be a matrix of logical, singles or doubles
    if (! all (all (Y == 0 | Y == 1)))
      error ("mnrfit: Y must contain only 1 and 0 when given as a 2D matrix.");
    endif
    ## Convert Y to a vector of positive integer categories
    Y = sum (bsxfun (@times, (1:K), Y), 2);
  endif

  ## Categorize Y in all other cases
  if (! iscellstr (Y))
    RowsUsed  = ! logical (sum (isnan ([X, Y]), 2));
    Y         = Y (RowsUsed);
    X         = X (RowsUsed, :);
    [UY, ~, YN] = unique (Y);  ## find unique categories in the response
    n = numel (UY);            ## number of unique response categories
  endif

  if (isnumeric (Y))
    if (! (all (Y > 0) && all (fix (Y) == Y)))
      error ("mnrfit: Y must contain positive integer category numbers.")
    endif
  endif

  ## Evaluate model type input argument
  switch (lower (MODELTYPE))
    case "nominal"
      if (n > 2)
        error ("mnrfit: fitting more than 2 nominal responses not supported.");
      else
        ## Y has two responses. Ordinal logistic regression can be used to fit
        ## models with binary nominal responses
      endif
    case "ordinal"
      ## Do nothing, ordinal responses are fully supported
    case "hierarchical"
      error ("mnrfit: fitting hierarchical responses not supported.");
    otherwise
      error ("mnrfit: model type not recognised.");
  endswitch

  ## Perform fit and reformat output
  [INTERCEPT, SLOPE, DEV, ~, ~, ~, S] = logistic_regression (YN - 1, X, dispopt);
  B = cat (1, INTERCEPT, SLOPE);
  STATS = struct ("beta", B, ...
                  "dfe", [], ...      ## Not used
                  "s", [], ...        ## Not used
                  "sfit", [], ...     ## Not used
                  "estdisp", [], ...  ## Not used
                  "coeffcorr", S.coeffcorr, ...
                  "covb", S.cov, ...
                  "se", S.se, ...
                  "t", [], ...        ## Not used
                  "p", [], ...        ## Not used
                  "resid", [], ...    ## Not used
                  "residp", [], ...   ## Not used
                  "residd", []);      ## Not used

endfunction

## Test input validation
%!error<mnrfit: too few input arguments.> mnrfit (ones (50,1))
%!error<mnrfit: Predictors must be numeric.> ...
%! mnrfit ({1 ;2 ;3 ;4 ;5}, ones (5,1))
%!error<mnrfit: Predictors must be a vector or a 2D matrix.> ...
%! mnrfit (ones (50, 4, 2), ones (50, 1))
%!error<mnrfit: Response must be a vector or a 2D matrix.> ...
%! mnrfit (ones (50, 4), ones (50, 1, 3))
%!error<mnrfit: Y must have the same number of rows as X.> ...
%! mnrfit (ones (50, 4), ones (45,1))
%!error<mnrfit: Response labels must be a character array, a cell vector> ...
%! mnrfit (ones (5, 4), {1 ;2 ;3 ;4 ;5})
%!error<mnrfit: optional arguments must be in pairs.> ...
%! mnrfit (ones (5, 4), ones (5, 1), "model")
%!error<mnrfit: Y must be a column vector when given as cellstr.> ...
%! mnrfit (ones (5, 4), {"q","q";"w","w";"q","q";"w","w";"q","q"})
%!error<mnrfit: Y must contain only 1 and 0 when given as a 2D matrix.> ...
%! mnrfit (ones (5, 4), [1, 2; 1, 2; 1, 2; 1, 2; 1, 2])
%!error<mnrfit: Y must contain positive integer category numbers.> ...
%! mnrfit (ones (5, 4), [1; -1; 1; 2; 1])
%!error<mnrfit: fitting more than 2 nominal responses not supported.> ...
%! mnrfit (ones (5, 4), [1; 2; 3; 2; 1], "model", "nominal")
%!error<mnrfit: fitting hierarchical responses not supported.> ...
%! mnrfit (ones (5, 4), [1; 2; 3; 2; 1], "model", "hierarchical")
%!error<mnrfit: model type not recognised.> ...
%! mnrfit (ones (5, 4), [1; 2; 3; 2; 1], "model", "whatever")

