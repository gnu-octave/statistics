## Copyright (C) 2020 Stefano Guidoni <ilguido@users.sf.net>
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{C} =} confusionmat (@var{group}, @var{grouphat})
## @deftypefnx {statistics} {@var{C} =} confusionmat (@var{group}, @var{grouphat}, "Order", @var{grouporder})
## @deftypefnx {statistics} {[@var{C}, @var{order}] =} confusionmat (@var{group}, @var{grouphat})
##
## Compute a confusion matrix for classification problems
##
## @code{confusionmat} returns the confusion matrix @var{C} for the group of
## actual values @var{group} and the group of predicted values @var{grouphat}.
## The row indices of the confusion matrix represent actual values, while the
## column indices represent predicted values. The indices are the same for both
## actual and predicted values, so the confusion matrix is a square matrix.
## Each element of the matrix represents the number of matches between a given
## actual value (row index) and a given predicted value (column index), hence
## correct matches lie on the main diagonal of the matrix.
## The order of the rows and columns is returned in @var{order}.
##
## @var{group} and @var{grouphat} must have the same number of observations
## and the same data type.
## Valid data types are numeric vectors, logical vectors, character arrays,
## string arrays (not implemented yet), cell arrays of strings.
##
## The order of the rows and columns can be specified by setting the
## @var{grouporder} variable. The data type of @var{grouporder} must be the
## same of @var{group} and @var{grouphat}.
##
## MATLAB compatibility: Octave misses string arrays and categorical vectors.
##
## @seealso{crosstab}
## @end deftypefn

function [C, order] = confusionmat (group, grouphat, opt = "Order", grouporder)

  ## check the input parameters
  if ((nargin < 2) || (nargin > 4))
    print_usage();
  endif

  y_true = group;
  y_pred = grouphat;

  if (class (y_true) != class (y_pred))
    error ("confusionmat: group and grouphat must be of the same data type");
  endif

  if (length (y_true) != length (y_pred))
    error ("confusionmat: group and grouphat must be of the same length");
  endif

  if ((nargin > 3) && strcmp (opt, "Order"))
    unique_tokens = grouporder;

    if class( y_true ) != class( unique_tokens )
      error ("confusionmat: group and grouporder must be of the same data type");
    endif
  endif

  if (isvector (y_true))
    if (isrow (y_true))
      y_true = vec(y_true);
    endif
  else
    error ("confusionmat: group must be a vector or array");
  endif

  if (isvector (y_pred))
    if (isrow (y_pred))
      y_pred = vec(y_pred);
    endif
  else
    error ("confusionmat: grouphat must be a vector or array");
  endif

  if (exist ( "unique_tokens", "var"))
    if (isvector (unique_tokens))
      if (isrow (unique_tokens))
        unique_tokens = vec(unique_tokens);
      endif
    else
      error ("confusionmat: grouporder must be a vector or array");
    endif
  endif

  ## compute the confusion matrix
  if (isa (y_true, "numeric") || isa (y_true, "logical"))
    ## numeric or boolean vector

    ## MATLAB compatibility:
    ## remove NaN values from grouphat
    nan_indices = find (isnan (y_pred));
    y_pred(nan_indices) = [];

    ## MATLAB compatibility:
    ## numeric and boolean values
    ## are sorted in ascending order
    if (! exist ("unique_tokens", "var"))
      unique_tokens = union (y_true, y_pred);
    endif

    y_true(nan_indices) = [];

    C_size = length (unique_tokens);

    C = zeros (C_size);

    for i = 1:length (y_true)
      row_index = find (unique_tokens == y_true(i));
      col_index = find (unique_tokens == y_pred(i));

      C(row_index, col_index)++;
    endfor

  elseif (iscellstr (y_true))
    ## string cells

    ## MATLAB compatibility:
    ## remove empty values from grouphat
    empty_indices = [];
    for i = 1:length (y_pred)
      if (isempty (y_pred{i}))
        empty_indices = [empty_indices; i];
      endif
    endfor

    y_pred(empty_indices) = [];

    ## MATLAB compatibility:
    ## string values are sorted according to their
    ## first appearance in group and grouphat
    if (! exist ("unique_tokens", "var"))
      all_tokens = vertcat (y_true, y_pred);
      unique_tokens = [all_tokens(1)];

      for i = 2:length( all_tokens )
        if (! any (strcmp (all_tokens(i), unique_tokens)))
          unique_tokens = [unique_tokens; all_tokens(i)];
        endif
      endfor
    endif

    y_true(empty_indices) = [];
    C_size = length (unique_tokens);
    C = zeros (C_size);

    for i = 1:length (y_true)
      row_index = find (strcmp (y_true{i}, unique_tokens));
      col_index = find (strcmp (y_pred{i}, unique_tokens));
      C(row_index, col_index)++;
    endfor

  elseif (ischar (y_true))
    ## character array

    ## MATLAB compatibility:
    ## character values are sorted according to their
    ## first appearance in group and grouphat
    if (! exist ("unique_tokens", "var"))
      all_tokens = vertcat (y_true, y_pred);
      unique_tokens = [all_tokens(1)];
      for i = 2:length (all_tokens)
        if (! any (find (unique_tokens == all_tokens(i))))
          unique_tokens = [unique_tokens; all_tokens(i)];
        endif
      endfor
    endif

    C_size = length ( unique_tokens );
    C = zeros ( C_size );

    for i = 1:length( y_true)
      row_index = find (unique_tokens == y_true(i));
      col_index = find (unique_tokens == y_pred(i));
      C(row_index, col_index)++;
    endfor

  // MY FIX STARTS
  elseif (isstring (y_true) || iscategorical (y_true))
    ## string or categorical array

    if (! exist ("unique_tokens", "var"))
      unique_tokens = union (y_true, y_pred);
    endif

    C_size = length (unique_tokens);
    C = zeros (C_size);

    for i = 1:length (y_true)
      row_index = find (unique_tokens == y_true(i));
      col_index = find (unique_tokens == y_pred(i));

      if (! isempty (row_index) && ! isempty (col_index))
        C(row_index, col_index)++;
      endif
    endfor
    // FIX COMPLETED
  else
    error ("confusionmat: invalid data type");
  endif

  order = unique_tokens;
endfunction

## Test the confusion matrix example from
## R.Bonnin, "Machine Learning for Developers", pp. 55-56
%!test
%! Yt = [8 5 6 8 5 3 1 6 4 2 5 3 1 4]';
%! Yp = [8 5 6 8 5 2 3 4 4 5 5 7 2 6]';
%! C  = [0 1 1 0 0 0 0 0; 0 0 0 0 1 0 0 0; 0 1 0 0 0 0 1 0; 0 0 0 1 0 1 0 0; ...
%!       0 0 0 0 3 0 0 0; 0 0 0 1 0 1 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 2];
%! assert (confusionmat (Yt, Yp), C)
