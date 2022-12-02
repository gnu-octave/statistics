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
## @deftypefn {Function File} @var{h} = levene_test (@var{x})
## @deftypefnx {Function File} @var{h} = levene_test (@var{x}, @var{group})
## @deftypefnx {Function File} @var{h} = levene_test (@var{x}, @var{alpha})
## @deftypefnx {Function File} @var{h} = levene_test (@var{x}, @var{testtype})
## @deftypefnx {Function File} @var{h} = levene_test (@var{x}, @var{group}, @var{alpha})
## @deftypefnx {Function File} @var{h} = levene_test (@var{x}, @var{group}, @var{testtype})
## @deftypefnx {Function File} @var{h} = levene_test (@var{x}, @var{group}, @var{alpha}, @var{testtype})
## @deftypefnx {Function File} [@var{h}, @var{pval}] = levene_test (@dots{})
## @deftypefnx {Function File} [@var{h}, @var{pval}, @var{W}] = levene_test (@dots{})
## @deftypefnx {Function File} [@var{h}, @var{pval}, @var{W}, @var{df}] = levene_test (@dots{})
##
## Perform a Levene's test for the homogeneity of variances.
##
## Under the null hypothesis of equal variances, the test statistic @var{W}
## approximately follows an F distribution with @var{df} degrees of
## freedom being a vector ([k-1, N-k]).
##
## The p-value (1 minus the CDF of this distribution at @var{W}) is returned in
## @var{pval}.  @var{h} = 1 if the null hypothesis is rejected at the
## significance level of @var{alpha}.  Otherwise @var{h} = 0.
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
##
## @item
## @var{testtype} is a string determining the type of Levene's test.  By default
## it is set to "absolute", but the user can also parse "quadratic" in order to
## perform Levene's Quadratic test for equal variances or "median" in order to
## to perform the Brown-Forsythe's test.  These options determine how the Z_ij
## values are computed.  If an invalid name is parsed for @var{testtype}, then
## the Levene's Absolute test is performed.
## @end itemize
##
## @seealso{bartlett_test, vartest2, vartestn}
## @end deftypefn

function [h, pval, W, df] = levene_test (x, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 4)
    error ("levene_test: invalid number of input arguments.");
  endif

  ## Add defaults
  group = [];
  alpha = 0.05;
  ttype = "absolute";

  ## Check for 2nd argument being ALPHA, GROUP, or TESTTYPE
  if (nargin > 1)
    if (isscalar (varargin{1}) && isnumeric (varargin{1}) ...
                               && numel (varargin{1}) == 1)
      alpha = varargin{1};
      ## Check for valid alpha value
      if (alpha <= 0 || alpha >= 1)
       error ("levene_test: wrong value for alpha.");
      endif
    elseif (any (strcmpi (varargin{1}, {"absolute", "quadratic", "median"})))
      ttype = varargin{1};
    elseif (isvector (varargin{1}) && numel (varargin{1} > 1))
      if ((size (x, 2) == 1 && size (x, 1) == numel (varargin{1})) || ...
          (size (x, 2) > 1 && size (x, 2) == numel (varargin{1})))
        group = varargin{1};
      else
        error ("levene_test: GROUP and X mismatch.");
      endif
    elseif (isempty (varargin{1}))
      ## Do nothing
    else
      error ("levene_test: invalid second input argument.");
    endif
  endif

  ## Check for 3rd argument
  if (nargin > 2)
    if (isscalar (varargin{2}) && isnumeric (varargin{2}) ...
                               && numel (varargin{2} == 1))
      alpha = varargin{2};
      ## Check for valid alpha value
      if (alpha <= 0 || alpha >= 1)
        error ("levene_test: wrong value for alpha.");
      endif
    elseif (any (strcmpi (varargin{2}, {"absolute", "quadratic", "median"})))
      ttype = varargin{2};
    else
      error ("levene_test: invalid third input argument.");
    endif
  endif

  ## Check for 3rd argument
  if (nargin > 3)
    if (any (strcmpi (varargin{3}, {"absolute", "quadratic", "median"})))
      ttype = varargin{3};
    else
      error ("levene_test: invalid option for TESTTYPE as 4th argument.");
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
      error ("levene_test: columns in X and GROUP length do not match.");
    endif
  endif

  ## Check that x and group are the same size
  if (! all (numel (x) == numel (group)))
    error (srtcat (["levene_test: GROUP must be a vector with the same"], ...
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

  ## Get sample size (N_i), mean (Y_i), median (Y_I) and sample values (Y_ij)
  ## for groups with more than one sample
  groups = size (group_names, 1);
  rgroup = [];
  N_i = zeros (1, groups);
  Y_i = N_i;
  Y_I = N_i;
  for k = 1:groups
    group_size = find (group_id == k);
    if (length (group_size) > 1)
      N_i(k) = length (group_size);
      Y_i(k) = mean (x(group_size));
      Y_I(k) = median (x(group_size));
      Y_ij{k} = x(group_size);
    else
      warning (strcat (sprintf ("levene_test: GROUP %s has a single", ...
              group_names{k}), [" sample and is not included in the test.\n"]));
      rgroup = [rgroup, k];
      N_i(k) = 1;
      Y_i(k) = x(group_size);
      Y_I(k) = x(group_size);
      Y_ij{k} = x(group_size);
    endif
  endfor

  ## Remove groups with a single sample
  if (! isempty (rgroup))
    N_i(rgroup) = [];
    Y_i(rgroup) = [];
    Y_I(rgroup) = [];
    Y_ij(rgroup) = [];
    k = k - numel (rgroup);
  endif

  ## Compute Z_ij for "absolute" or "quadratic" Levene's test
  switch (lower (ttype))
    case "absolute"
      for i = 1:k
        Z_ij{i} = abs (Y_ij{i} - Y_i(i));
      endfor
    case "quadratic"
      for i = 1:k
        Z_ij{i} = sqrt ((Y_ij{i} - Y_i(i)) .^ 2);
      endfor
    case "median"
      for i = 1:k
        Z_ij{i} = abs (Y_ij{i} - Y_I(i));
      endfor
  endswitch

  ## Compute Z_i and Z_
  Z_ = [];
  for i = 1:k
    Z_i(i) = mean (Z_ij{i});
    Z_ = [Z_; Z_ij{i}(:)];
  endfor
  Z_ = mean (Z_);

  ## Compute total sample size (N)
  N = sum (N_i);

  ## Calculate W statistic.
  termA = (N - k) / (k - 1);
  termB = sum (N_i .* ((Z_i - Z_) .^ 2));
  termC = 0;
  for i = 1:k
    termC += sum ((Z_ij{i} - Z_i(i)) .^ 2);
  endfor
  W = termA * (termB / termC);

  ## Calculate p-value from the chi-square distribution
  pval = 1 - fcdf (W, k - 1, N - k);

  ## Save dfs
  df = [k-1, N-k];

  ## Determine the test outcome
  h = double (pval < alpha);
endfunction

## Test input validation
%!error<levene_test: invalid number of input arguments.> levene_test ()
%!error<levene_test: invalid number of input arguments.> ...
%! levene_test (1, 2, 3, 4, 5);
%!error<levene_test: wrong value for alpha.> levene_test (randn (50, 2), 0);
%!error<levene_test: GROUP and X mismatch.> ...
%! levene_test (randn (50, 2), [1, 2, 3]);
%!error<levene_test: GROUP and X mismatch.> ...
%! levene_test (randn (50, 1), ones (55, 1));
%!error<levene_test: invalid second input argument.> ...
%! levene_test (randn (50, 1), ones (50, 2));
%!error<levene_test: wrong value for alpha.> ...
%! levene_test (randn (50, 2), [], 1.2);
%!error<levene_test: GROUP and X mismatch.> ...
%! levene_test (randn (50, 2), "some_string");
%!error<levene_test: invalid third input argument.> ...
%! levene_test (randn (50, 2), [], "alpha");
%!error<levene_test: wrong value for alpha.> ...
%! levene_test (randn (50, 1), [ones(25, 1); 2*ones(25, 1)], 1.2);
%!error<levene_test: invalid third input argument.> ...
%! levene_test (randn (50, 1), [ones(25, 1); 2*ones(25, 1)], "err");
%!error<levene_test: invalid option for TESTTYPE as 4th argument.> ...
%! levene_test (randn (50, 1), [ones(25, 1); 2*ones(25, 1)], 0.05, "type");
%!warning<levene_test: GROUP> ...
%! levene_test (randn (50, 1), [ones(24, 1); 2*ones(25, 1); 3]);
## Test results
%!test
%! load examgrades
%! [h, pval, W, df] = levene_test (grades);
%! assert (h, 1);
%! assert (pval, 9.523239714592791e-07, 1e-14);
%! assert (W, 8.59529, 1e-5);
%! assert (df, [4, 595]);
%!test
%! load examgrades
%! [h, pval, W, df] = levene_test (grades, [], "quadratic");
%! assert (h, 1);
%! assert (pval, 9.523239714592791e-07, 1e-14);
%! assert (W, 8.59529, 1e-5);
%! assert (df, [4, 595]);
%!test
%! load examgrades
%! [h, pval, W, df] = levene_test (grades, [], "median");
%! assert (h, 1);
%! assert (pval, 1.312093241723211e-06, 1e-14);
%! assert (W, 8.415969, 1e-6);
%! assert (df, [4, 595]);
%!test
%! load examgrades
%! [h, pval, W, df] = levene_test (grades(:,[1:3]));
%! assert (h, 1);
%! assert (pval, 0.004349390980463497, 1e-14);
%! assert (W, 5.52139, 1e-5);
%! assert (df, [2, 357]);
%!test
%! load examgrades
%! [h, pval, W, df] = levene_test (grades(:,[1:3]), "median");
%! assert (h, 1);
%! assert (pval, 0.004355216763951453, 1e-14);
%! assert (W, 5.52001, 1e-5);
%! assert (df, [2, 357]);
%!test
%! load examgrades
%! [h, pval, W, df] = levene_test (grades(:,[3,4]), "quadratic");
%! assert (h, 0);
%! assert (pval, 0.1807494957440653, 2e-14);
%! assert (W, 1.80200, 1e-5);
%! assert (df, [1, 238]);
%!test
%! load examgrades
%! [h, pval, W, df] = levene_test (grades(:,[3,4]), "median");
%! assert (h, 0);
%! assert (pval, 0.1978225622063785, 2e-14);
%! assert (W, 1.66768, 1e-5);
%! assert (df, [1, 238]);
