## Copyright (C) 1995-2017 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {[@var{h}, @var{pval}, @var{stats}] =} regression_ftest (@var{y}, @var{x}, @var{fm})
## @deftypefnx {statistics} {[@dots{}] =} regression_ftest (@var{y}, @var{x}, @var{fm}, @var{rm})
## @deftypefnx {statistics} {[@dots{}] =} regression_ftest (@var{y}, @var{x}, @var{fm}, @var{rm}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@dots{}] =} regression_ftest (@var{y}, @var{x}, @var{fm}, [], @var{Name}, @var{Value})
##
## F-test for General Linear Regression Analysis
##
## Perform a general linear regression F test for the null hypothesis that the
## full model of the form @nospell{y = b_0 + b_1 * x_1 + b_2 * x_2 + @dots{} +
## b_n * x_n + e}, where n is the number of variables in @var{x}, does not
## perform better than a reduced model, such as @nospell{y = b'_0 + b'_1 * x_1 +
## b'_2 * x_2 + @dots{} + b'_k * x_k + e}, where k < n and it corresponds to the
## first k variables in @var{x}.  Explanatory (dependent) variable @var{y} and
## response (independent) variables @var{x} must not contain any missing values
## (NaNs).
##
## The full model, @var{fm}, must be a vector of length equal to the columns of
## @var{x}, in which case the constant term b_0 is assumed 0, or equal to
## the columns of @var{x} plus one, in which case the first element is the
## constant b_0.
##
## The reduced model, @var{rm}, must include the constant term and a subset of
## the variables (columns) in @var{x}. If @var{rm} is not given, then a constant
## term b'_0 is assumed equal to the constant term, b_0, of the full model or 0,
## if the full model, @var{fm}, does not have a constant term.  @var{rm} must be
## a vector or a scalar if only a constant term is passed into the function.
##
## Name-Value pair arguments can be used to set statistical significance.
## @qcode{"alpha"} can be used to specify the significance level of the test
## (the default value is 0.05).  If you want pass optional Name-Value pair
## without a reduced model, make sure that the latter is passed as an empty
## variable.
##
## If @var{h} is 1 the null hypothesis is rejected, meaning that the full model
## explains the variance better than the restricted model.  If @var{h} is 0, it
## can be assumed that the full model does NOT explain the variance any better
## than the restricted model.
##
## The p-value (1 minus the CDF of this distribution at @var{f}) is returned
## in @var{pval}.
##
## Under the null, the test statistic @var{f} follows an F distribution with
## 'df1' and 'df2' degrees of freedom, which are returned as fields in the
## @var{stats} structure along with the test's F-statistic, 'fstat'
##
## @seealso{regress, regression_ttest}
## @end deftypefn

function [h, pval, stats] = regression_ftest (y, x, fm, rm, varargin)

  ## Check for valid input
  if (nargin < 3)
    print_usage ();
  endif

  ## Check for finite real numbers in Y, X
  if (! all (isfinite (y)) || ! isreal (y))
    error ("regression_ftest: Y must contain finite real numbers.");
  endif
  if (! all (isfinite (x(:))) || ! isreal (x))
    error ("regression_ftest: X must contain finite real numbers.");
  endif

  ## Set default arguments
  alpha = 0.05;

  ## Check additional options
  i = 1;
  while (i <= length (varargin))
    switch lower (varargin{i})
      case "alpha"
        i = i + 1;
        alpha = varargin{i};
        ## Check for valid alpha
        if (! isscalar (alpha) || ! isnumeric (alpha) || ...
                    alpha <= 0 || alpha >= 1)
          error ("regression_ftest: invalid value for alpha.");
        endif
      otherwise
        error ("regression_ftest: invalid Name argument.");
    endswitch
    i = i + 1;
  endwhile

  ## Get size of response (independent) variables
  [s, v] = size (x);
  ## Add a constant term of 1s in X
  x = [ones(s, 1), x];

  ## Check the size of explanatory (dependent) variable
  if (! (isvector (y) && (length (y) == s)))
    error ("regression_ftest: Y must be a vector of length 'rows (X)'.");
  endif
  y = reshape (y, s, 1);

  ## Check the full model
  if (! (isvector (fm) && (length (fm) == v || length (fm) == v + 1)))
    error (strcat (["regression_ftest: full model, FM, must be a vector"], ...
                   [" of length equal to 'rows (X)' or 'rows (X) + 1'."]));
  endif
  ## Make it row vector and add a constant = 0 if necessary
  fm_len = length (fm);
  fm = reshape (fm, 1, fm_len);
  if (fm_len == v)
    fm = [0, fm];
    fm_len += 1;
  endif

  ## Check the reduced model
  if (nargin - length (varargin) == 4)
    if (isempty (rm))
      rm = [fm(1), zeros(1, fm_len - 1)];
      rm_len = 1;
    else
      if (! isvector (rm) || ! isnumeric (rm))
        error (strcat (["regression_ftest: reduced model, RM, must be a"], ...
                       [" numeric vector or a scalar."]));
      endif
      rm_len = length (rm);
      if (rm_len >= fm_len - 1)
        error (strcat (["regression_ftest: reduced model, RM, must have"], ...
                       [" smaller length than the full model, FM."]));
      endif
      rm = reshape (rm, 1, rm_len);
      rm = [rm, zeros(1, fm_len - rm_len)];
    endif
  else
    rm = [fm(1), zeros(1, fm_len - 1)];
    rm_len = 1;
  endif

  ## Calculate the fitted response for full and reduced models
  y_fm = sum (x .* fm, 2);
  y_rm = sum (x .* rm, 2);

  ## Calculate Sum of Squares Error for full and reduced models
  SSE_fm = sumsq (y - y_fm);
  SSE_rm = sumsq (y - y_rm);

  ## Calculate the necessary statistics
  stats.df1 = fm_len - rm_len;
  stats.df2 = s - v;
  stats.fstat = ((SSE_rm - SSE_fm) / stats.df1) / (SSE_fm / stats.df2);
  pval = 1 - fcdf (stats.fstat, stats.df1, stats.df2);

  ## Determine the test outcome
  ## MATLAB returns this a double instead of a logical array
  h = double (pval < alpha);

endfunction

## Test input validation
%!error<Invalid call to regression_ftest.  Correct usage> regression_ftest ();
%!error<Invalid call to regression_ftest.  Correct usage> ...
%! regression_ftest ([1 2 3]', [2 3 4; 3 4 5]');
%!error<regression_ftest: Y must contain finite real numbers.> ...
%! regression_ftest ([1 2 NaN]', [2 3 4; 3 4 5]', [1 0.5]);
%!error<regression_ftest: Y must contain finite real numbers.> ...
%! regression_ftest ([1 2 Inf]', [2 3 4; 3 4 5]', [1 0.5]);
%!error<regression_ftest: Y must contain finite real numbers.> ...
%! regression_ftest ([1 2 3+i]', [2 3 4; 3 4 5]', [1 0.5]);
%!error<regression_ftest: X must contain finite real numbers.> ...
%! regression_ftest ([1 2 3]', [2 3 NaN; 3 4 5]', [1 0.5]);
%!error<regression_ftest: X must contain finite real numbers.> ...
%! regression_ftest ([1 2 3]', [2 3 Inf; 3 4 5]', [1 0.5]);
%!error<regression_ftest: X must contain finite real numbers.> ...
%! regression_ftest ([1 2 3]', [2 3 4; 3 4 3+i]', [1 0.5]);
%!error<regression_ftest: invalid value for alpha.> ...
%! regression_ftest ([1 2 3]', [2 3 4; 3 4 5]', [1 0.5], [], "alpha", 0);
%!error<regression_ftest: invalid value for alpha.> ...
%! regression_ftest ([1 2 3]', [2 3 4; 3 4 5]', [1 0.5], [], "alpha", 1.2);
%!error<regression_ftest: invalid value for alpha.> ...
%! regression_ftest ([1 2 3]', [2 3 4; 3 4 5]', [1 0.5], [], "alpha", [.02 .1]);
%!error<regression_ftest: invalid value for alpha.> ...
%! regression_ftest ([1 2 3]', [2 3 4; 3 4 5]', [1 0.5], [], "alpha", "a");
%!error<regression_ftest: invalid Name argument.> ...
%! regression_ftest ([1 2 3]', [2 3 4; 3 4 5]', [1 0.5], [], "some", 0.05);
%!error<regression_ftest: Y must be a vector of length> ...
%! regression_ftest ([1 2 3]', [2 3; 3 4]', [1 0.5]);
%!error<regression_ftest: Y must be a vector of length> ...
%! regression_ftest ([1 2; 3 4]', [2 3; 3 4]', [1 0.5]);
%!error<regression_ftest: reduced model, RM, must be a numeric vector or> ...
%! regression_ftest ([1 2 3]', [2 3 4; 3 4 5]', [1 0.5], ones (2));
%!error<regression_ftest: reduced model, RM, must be a numeric vector or> ...
%! regression_ftest ([1 2 3]', [2 3 4; 3 4 5]', [1 0.5], "alpha");
%!error<regression_ftest: reduced model, RM, must have smaller length than> ...
%! regression_ftest ([1 2 3]', [2 3 4; 3 4 5]', [1 0.5], [1 2]);

## Test results
