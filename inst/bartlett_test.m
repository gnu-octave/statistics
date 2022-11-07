## Copyright (C) 1995-2017 Kurt Hornik
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} @var{h} = bartlett_test (@var{x})
## @deftypefnx {Function File} @var{h} = bartlett_test (@var{x}, @var{group})
## @deftypefnx {Function File} @var{h} = bartlett_test (@var{x}, @var{alpha})
## @deftypefnx {Function File} @var{h} = bartlett_test (@var{x}, @var{group}, @var{alpha})
## @deftypefnx {Function File} [@var{h}, @var{pval}] = bartlett_test (@dots{})
## @deftypefnx {Function File} [@var{h}, @var{pval}, @var{chisq}] = bartlett_test (@dots{})
## @deftypefnx {Function File} [@var{h}, @var{pval}, @var{chisq}, @var{df}] = bartlett_test (@dots{})
##
## Perform a Bartlett test for the homogeneity of variances.
##
## Under the null hypothesis of equal variances, the test statistic @var{chisq}
## approximately follows a chi-square distribution with @var{df} degrees of
## freedom.
##
## The p-value (1 minus the CDF of this distribution at @var{chisq}) is
## returned in @var{pval}.  @var{h} = 1 if the null hypothesis is rejected at
## the significance level of @var{alpha}.  Otherwise @var{h} = 0.
##
## Input Arguments:
##
## @itemize
## @item
## @var{x} contains the data and it can either be a vector or matrix.
## If @var{x} is a matrix, then each column is treated as a separate group.
## If @var{x} is a vector, then the @var{group} argument is mandatory.
## NaN values are omitted.
##
## @item
## @var{group} contains the names for each group.  If @var{x} is a vector, then
## @var{group} must be a vector of the same length, or a string array or cell
## array of strings with one row for each element of @var{x}.  @var{x} values
## corresponding to the same value of @var{group} are placed in the same group.
## If @var{x} is a matrix, then @var{group} can either be a cell array of
## strings of a character array, with one row per column of @var{x} in the same
## way it is used in @code{anova1} function.  If @var{x} is a matrix, then
## @var{group} can be omitted either by entering an empty array ([]) or by
## parsing only @var{alpha} as a second argument (if required to change its
## default value).
##
## @item
## @var{alpha} is the statistical significance value at which the null
## hypothesis is rejected.  Its default value is 0.05 and it can be parsed
## either as a second argument (when @var{group} is omitted) or as a third
## argument.
## @end itemize
##
## @seealso{levene_test, vartest2, vartestn}
## @end deftypefn

function [h, pval, chisq, df] = bartlett_test (x, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    error ("bartlett_test: invalid number of input arguments.");
  endif

  ## Add defaults
  group = [];
  alpha = 0.05;

  ## Check for 2nd argument being ALPHA or GROUP
  if (nargin > 1)
    if (isscalar (varargin{1}) && isnumeric (varargin{1}) ...
                               && numel (varargin{1}) == 1)
      alpha = varargin{1};
      ## Check for valid alpha value
      if (alpha <= 0 || alpha >= 1)
       error ("bartlett_test: wrong value for alpha.");
      endif
    elseif (isvector (varargin{1}) && numel (varargin{1} > 1))
      if ((size (x, 2) == 1 && size (x, 1) == numel (varargin{1})) || ...
          (size (x, 2) > 1 && size (x, 2) == numel (varargin{1})))
        group = varargin{1};
      else
        error ("bartlett_test: GROUP and X mismatch.");
      endif
    elseif (isempty (varargin{1}))
      ## Do nothing
    else
      error ("bartlett_test: invalid second input argument.");
    endif
  endif

  ## Check for 3rd argument
  if (nargin > 2)
    alpha = varargin{2};
    ## Check for valid alpha value
    if (! isscalar (alpha) || ! isnumeric (alpha) || alpha <= 0 || alpha >= 1)
      error ("bartlett_test: wrong value for alpha.");
    endif
  endif

  ## Convert group to cell array from character array, make it a column
  if (! isempty (group) && ischar (group))
    group = cellstr (group);
  endif
  if (size (group, 1) == 1)
    group = group';
  endif

  ## If x is a matrix, convert it to column vector and create a
  ## corresponging column vector for groups
  if (length (x) < prod (size (x)))
    [n, m] = size (x);
    x = x(:);
    gi = reshape (repmat ((1:m), n, 1), n*m, 1);
    if (length (group) == 0)          ## no group names are provided
      group = gi;
    elseif (size (group, 1) == m)     ## group names exist and match columns
      group = group(gi,:);
    else
      error ("bartlett_test: columns in X and GROUP length do not match.");
    endif
  endif

  ## Check that x and group are the same size
  if (! all (numel (x) == numel (group)))
    error (srtcat (["bartlett_test: GROUP must be a vector with the same"], ...
                   [" number of rows as x."]));
  endif

  ## Identify NaN values (if any) and remove them from X along with
  ## their corresponding values from group vector
  nonan = ! isnan (x);
  x = x(nonan);
  group = group(nonan, :);

  ## Convert group to indices and separate names
  [group_id, group_names] = grp2idx (group);
  group_id = group_id(:);

  ## Get sample size (n_i) and var (s^2_i) for each group with n_i > 1
  groups = size (group_names, 1);
  rgroup = [];
  n_i = zeros (1, groups);
  s_i = n_i;
  for k = 1:groups
    group_size = find (group_id == k);
    if (length (group_size) > 1)
      n_i(k) = length (group_size);
      s_i(k) = var (x(group_size));
    else
      warning (strcat (sprintf ("bartlett_test: GROUP %s has a single", ...
              group_names{k}), [" sample and is not included in the test.\n"]));
      rgroup = [rgroup, k];
      n_i(k) = 1;
      s_i(k) = NaN;
    endif
  endfor

  ## Remove groups with a single sample
  if (! isempty (rgroup))
    n_i(rgroup) = [];
    s_i(rgroup) = [];
    k = k - numel (rgroup);
  endif

  ## Compute total sample size (N) and pooled variance (S)
  N = sum (n_i);
  S = (1 / (N - k)) * sum ((n_i - 1) .* s_i);

  ## Calculate B statistic. That is, B ~ X^2(k-1)
  B_nom = (N - k) * log (S) - sum ((n_i - 1) .* log (s_i));
  B_den = 1 + (1 / (3 * (k - 1))) * (sum (1 ./ (n_i - 1)) - (1 / (N - k)));
  chisq = B_nom / B_den;

  ## Calculate p-value from the chi-square distribution
  df = k - 1;
  pval = 1 - chi2cdf (chisq, df);

  ## Determine the test outcome
  h = double (pval < alpha);
endfunction

## Test input validation
%!error<bartlett_test: invalid number of input arguments.> bartlett_test ()
%!error<bartlett_test: invalid number of input arguments.> ...
%! bartlett_test (1, 2, 3, 4);
%!error<bartlett_test: wrong value for alpha.> bartlett_test (randn (50, 2), 0);
%!error<bartlett_test: GROUP and X mismatch.> ...
%! bartlett_test (randn (50, 2), [1, 2, 3]);
%!error<bartlett_test: GROUP and X mismatch.> ...
%! bartlett_test (randn (50, 1), ones (55, 1));
%!error<bartlett_test: invalid second input argument.> ...
%! bartlett_test (randn (50, 1), ones (50, 2));
%!error<bartlett_test: wrong value for alpha.> ...
%! bartlett_test (randn (50, 2), [], 1.2);
%!error<bartlett_test: wrong value for alpha.> ...
%! bartlett_test (randn (50, 2), [], "alpha");
%!error<bartlett_test: wrong value for alpha.> ...
%! bartlett_test (randn (50, 1), [ones(25, 1); 2*ones(25, 1)], 1.2);
%!error<bartlett_test: wrong value for alpha.> ...
%! bartlett_test (randn (50, 1), [ones(25, 1); 2*ones(25, 1)], "err");
%!warning<bartlett_test: GROUP> ...
%! bartlett_test (randn (50, 1), [ones(24, 1); 2*ones(25, 1); 3]);
## Test results
%!test
%! load examgrades
%! [h, pval, chisq, df] = bartlett_test (grades);
%! assert (h, 1);
%! assert (pval, 7.908647337018238e-08, 1e-14);
%! assert (chisq, 38.73324, 1e-5);
%! assert (df, 4);
%!test
%! load examgrades
%! [h, pval, chisq, df] = bartlett_test (grades(:,[2:4]));
%! assert (h, 1);
%! assert (pval, 0.01172, 1e-5);
%! assert (chisq, 8.89274, 1e-5);
%! assert (df, 2);
%!test
%! load examgrades
%! [h, pval, chisq, df] = bartlett_test (grades(:,[1,4]));
%! assert (h, 0);
%! assert (pval, 0.88118, 1e-5);
%! assert (chisq, 0.02234, 1e-5);
%! assert (df, 1);
%!test
%! load examgrades
%! grades = [grades; nan(10, 5)];
%! [h, pval, chisq, df] = bartlett_test (grades(:,[1,4]));
%! assert (h, 0);
%! assert (pval, 0.88118, 1e-5);
%! assert (chisq, 0.02234, 1e-5);
%! assert (df, 1);
%!test
%! load examgrades
%! [h, pval, chisq, df] = bartlett_test (grades(:,[2,5]), 0.01);
%! assert (h, 0);
%! assert (pval, 0.01791, 1e-5);
%! assert (chisq, 5.60486, 1e-5);
%! assert (df, 1);
