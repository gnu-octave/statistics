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
## string arrays, cell arrays of strings, and categorical arrays.
##
## The order of the rows and columns can be specified by setting the
## @var{grouporder} variable. The data type of @var{grouporder} must be the
## same of @var{group} and @var{grouphat}.
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

  obs_true = length (y_true);
  if (ischar (y_true))
    obs_true = rows (y_true);
  endif

  obs_pred = length (y_pred);
  if (ischar (y_pred))
    obs_pred = rows (y_pred);
  endif

  if (obs_true != obs_pred)
    error ("confusionmat: group and grouphat must be of the same length");
  endif

  if ((nargin > 3) && strcmp (opt, "Order"))
    unique_tokens = grouporder;

    if class( y_true ) != class( unique_tokens )
      error ("confusionmat: group and grouporder must be of the same data type");
    endif
  endif

  convert_order_to_char = false;

  ## Handle y_true
  if (ischar (y_true))
    y_true = cellstr (y_true);
    convert_order_to_char = true;
  elseif (isvector (y_true))
    y_true = y_true(:);
  else
    error ("confusionmat: group must be a vector or character array");
  endif

  ## Handle y_pred
  if (ischar (y_pred))
    y_pred = cellstr (y_pred);
  elseif (isvector (y_pred))
    y_pred = y_pred(:);
  else
    error ("confusionmat: grouphat must be a vector or character array");
  endif

  if (exist ( "unique_tokens", "var"))
    if (ischar (unique_tokens))
      unique_tokens = cellstr (unique_tokens);
    elseif (isvector (unique_tokens))
      unique_tokens = unique_tokens(:);
    else
      error ("confusionmat: grouporder must be a vector or character array");
    endif
  endif

  ## compute the confusion matrix
  if (isa (y_true, "numeric") || isa (y_true, "logical"))

    ## Remove observations where EITHER group or grouphat is NaN
    nan_indices = isnan (y_true) | isnan (y_pred);
    y_true(nan_indices) = [];
    y_pred(nan_indices) = [];

    ## numeric and boolean values are sorted in ascending order
    if (! exist ("unique_tokens", "var"))
      unique_tokens = unique ([y_true; y_pred]);
    endif

    C_size = length (unique_tokens);
    C = zeros (C_size);

    for i = 1:length (y_true)
      row_index = find (unique_tokens == y_true(i));
      col_index = find (unique_tokens == y_pred(i));

      ## Check valid indices
      if (!isempty(row_index) && !isempty(col_index))
        C(row_index, col_index)++;
      endif
    endfor

  elseif (iscellstr (y_true))
    ## string cells

    ## remove observations where EITHER input is empty
    empty_indices = cellfun ("isempty", y_true) | cellfun ("isempty", y_pred);
    y_true(empty_indices) = [];
    y_pred(empty_indices) = [];

    ## string values are sorted according to their
    ## first appearance in group and grouphat
    if (! exist ("unique_tokens", "var"))
      all_tokens = [y_true; y_pred];
      if (isempty (all_tokens))
        unique_tokens = {};
      else
        [~, idx] = unique (all_tokens, "first");
        unique_tokens = all_tokens(sort (idx));
      endif
    endif

    C_size = length (unique_tokens);
    C = zeros (C_size);

    for i = 1:length (y_true)
      row_index = find (strcmp (y_true{i}, unique_tokens));
      col_index = find (strcmp (y_pred{i}, unique_tokens));
      
      if (!isempty(row_index) && !isempty(col_index))
        C(row_index, col_index)++;
      endif
    endfor

  elseif (ischar (y_true))

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

  elseif (isa (y_true, "string"))

    ## 1. Filter Missing Values
    bad_indices = ismissing (y_true) | ismissing (y_pred);
    
    y_true(bad_indices) = [];
    y_pred(bad_indices) = [];

    ## 2. Determine Order
    ## String arrays are sorted ALPHABETICALLY by unique().
    if (! exist ("unique_tokens", "var"))
      all_tokens = [y_true; y_pred];
      
      if (isempty (all_tokens))
         unique_tokens = strings (0, 1);
      else
         unique_tokens = unique (all_tokens);
      endif
    endif

    C_size = length (unique_tokens);
    C = zeros (C_size);

    [~, row_indices] = ismember (y_true, unique_tokens);
    [~, col_indices] = ismember (y_pred, unique_tokens);

    valid_mask = (row_indices > 0) & (col_indices > 0);
    row_indices = row_indices(valid_mask);
    col_indices = col_indices(valid_mask);

    for i = 1:length (row_indices)
      C(row_indices(i), col_indices(i))++;
    endfor

  elseif (isa (y_true, "categorical"))

    ## 1. Filter Undefined Values
    bad_indices = isundefined (y_true) | isundefined (y_pred);
    
    y_true(bad_indices) = [];
    y_pred(bad_indices) = [];

    ## 2. Determine Order
    if (! exist ("unique_tokens", "var"))
       ## This ensures the matrix includes all defined categories 
       cats_true = categories (y_true);
       cats_pred = categories (y_pred);
       
       ## Union of defined categories 
       all_cats = union (cats_true, cats_pred, "stable");
       
       ## Create the reference order vector
       unique_tokens = categorical (all_cats, all_cats);
    endif

    C_size = length (unique_tokens);
    C = zeros (C_size);

    [~, row_indices] = ismember (y_true, unique_tokens);
    [~, col_indices] = ismember (y_pred, unique_tokens);
    
    valid_mask = (row_indices > 0) & (col_indices > 0);
    row_indices = row_indices(valid_mask);
    col_indices = col_indices(valid_mask);

    for i = 1:length (row_indices)
      C(row_indices(i), col_indices(i))++;
    endfor

  else
    error ("confusionmat: invalid data type");
  endif

  order = unique_tokens;

  if (convert_order_to_char && iscellstr (order))
    order = char (order);
  endif

endfunction

## Test 1: Standard Numeric Vector (R.Bonnin example)
%!test
%! Yt = [8 5 6 8 5 3 1 6 4 2 5 3 1 4]';
%! Yp = [8 5 6 8 5 2 3 4 4 5 5 7 2 6]';
%! C  = [0 1 1 0 0 0 0 0; 0 0 0 0 1 0 0 0; 0 1 0 0 0 0 1 0; 0 0 0 1 0 1 0 0; ...
%!       0 0 0 0 3 0 0 0; 0 0 0 1 0 1 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 2];
%! assert (confusionmat (Yt, Yp), C)
## Test 2: Basic Integers
%!test
%! g  = [1; 2; 3; 1];
%! gh = [1; 2; 2; 1];
%! [C, order] = confusionmat (g, gh);
%! assert (C, [2 0 0; 0 1 0; 0 1 0]);
%! assert (order, [1; 2; 3]);
## Test 3: Logical Vectors
%!test
%! g  = [true; false; true; false];
%! gh = [true; true;  false; false];
%! [C, order] = confusionmat (g, gh);
%! assert (C, [1 1; 1 1]);
%! assert (order, [false; true]);
## Test 4: Floating Point Numbers
%!test
%! g  = [1.1; 2.2; 1.1];
%! gh = [1.1; 2.2; 2.2];
%! [C, order] = confusionmat (g, gh);
%! assert (C, [1 1; 0 1]);
%! assert (order, [1.1; 2.2]);
## Test 5: Numeric with NaNs 
%!test
%! g  = [1; 2; NaN; 3];
%! gh = [1; 1; 2;   3];
%! [C, order] = confusionmat (g, gh);
%! assert (C, [1 0 0; 1 0 0; 0 0 1]);
%! assert (order, [1; 2; 3]);
## Test 6: Empty Inputs
%!error
%! confusionmat ([], [])
## Test 7: Scalar Inputs 
%!test
%! [C, order] = confusionmat (1, 1);
%! assert (C, 1);
%! assert (order, 1);
## Test 8: Cell Array with Empty Strings 
%!test
%! g  = {'A'; ''; 'B'};
%! gh = {'A'; 'B'; 'B'};
%! [C, order] = confusionmat (g, gh);
%! assert (C, [1 0; 0 1]);
%! assert (order, {'A'; 'B'});
## Test 9: Character Arrays
%!test
%! g  = ['AA'; 'BB'; 'AA'; 'CC'];
%! gh = ['AA'; 'BB'; 'BB'; 'CC'];
%! [C, order] = confusionmat (g, gh);
%! assert (C, [1 1 0; 0 1 0; 0 0 1]);
%! assert (order, ['AA'; 'BB'; 'CC']);
## Test 10: Character Arrays (Whitespace Handling)
%!test
%! g  = char ('A', 'B', 'A');
%! gh = char ('A', 'A', 'B');
%! [C, order] = confusionmat (g, gh);
%! assert (C, [1 1; 1 0]);
%! assert (order, char ('A', 'B'));
## Test 11: Cell Array of Strings
%!test
%! g  = {'Cat'; 'Dog'; 'Cat'; 'Bird'};
%! gh = {'Cat'; 'Cat'; 'Bird'; 'Bird'};
%! [C, order] = confusionmat (g, gh);
%! assert (C, [1 0 1; 1 0 0; 0 0 1]);
%! assert (order, {'Cat'; 'Dog'; 'Bird'});
## Test 12: String Arrays 
%!test
%! g  = ["Apple"; "Banana"; "Apple"];
%! gh = ["Apple"; "Apple";  "Cherry"];
%! [C, order] = confusionmat (g, gh);
%! assert (C, [1 0 1; 1 0 0; 0 0 0]);
%! assert (order, ["Apple"; "Banana"; "Cherry"]);
## Test 13: String Arrays (Missing Values)
%!test
%! g = string ({"A"; "B"; "B"});
%! g(2) = missing;
%! gh = string(["A"; "B"; "B"]);
%! [C, order] = confusionmat (g, gh);
%! assert (C, [1 0; 0 1]);
%! assert (isequal (order, string(["A"; "B"])));
## Test 14: Categorical Arrays
%!test
%! g  = categorical ({'Small', 'Medium', 'Large'});
%! gh = categorical ({'Small', 'Large',  'Large'});
%! [C, order] = confusionmat (g, gh);
%! assert (C, [1 0 0; 1 0 0; 0 0 1]);
%! assert (cellstr (char (order)), {'Large'; 'Medium'; 'Small'});
## Test 15: Categorical (Undefined Values / NaN)
%!test
%! g = categorical ({'Red', 'Blue', 'Red'});
%! g(2) = NaN;
%! gh = categorical ({'Red', 'Blue', 'Red'});
%! [C, order] = confusionmat (g, gh);
%! assert (C, [0 0; 0 2]);
%! assert (cellstr (char (order)), {'Blue'; 'Red'});
## Test 16: Categorical (Unused Categories)
%!test
%! vals = {'A', 'B', 'A'};
%! cats = {'A', 'B', 'C'};
%! g  = categorical (vals, cats);
%! gh = categorical (vals, cats);
%! [C, order] = confusionmat (g, gh);
%! assert (size (C), [3 3]);
%! assert (C(3,3), 0);
%! assert (cellstr (char (order)), {'A'; 'B'; 'C'});
## Test 17: Categorical (Union of Categories)
%!test
%! g  = categorical ({'A'}, {'A', 'B'});
%! gh = categorical ({'A'}, {'A', 'C'});
%! [C, order] = confusionmat (g, gh);
%! assert (size (C), [3 3]);
%! assert (cellstr (char (order)), {'A'; 'B'; 'C'});
## Test 18: Row vs Column Vector
%!test
%! g  = [1, 2, 3];
%! gh = [1; 2; 3];
%! [C, order] = confusionmat (g, gh);
%! assert (C, eye(3));
%! assert (order, [1; 2; 3]);
## Test 19: Custom Order 
%!test
%! g  = [1; 2; 3];
%! gh = [1; 2; 3];
%! myOrder = [3; 2; 1];
%! [C, order] = confusionmat (g, gh, "Order", myOrder);
%! assert (C, [1 0 0; 0 1 0; 0 0 1]);
%! assert (order, [3; 2; 1]);
## Test 20: Custom Order (Reordering Strings)
%!test
%! g  = {'A'; 'B'};
%! gh = {'A'; 'B'};
%! [C, order] = confusionmat (g, gh, "Order", {'B'; 'A'});
%! assert (C, [1 0; 0 1]);
%! assert (order, {'B'; 'A'});
## Test 21: Custom Order (Subset / Filtering)
%!test
%! g  = [1; 2; 3];
%! gh = [1; 2; 3];
%! [C, order] = confusionmat (g, gh, "Order", [1; 2]);
%! assert (C, eye(2));
%! assert (order, [1; 2]);
## Test 22: Custom Order (Superset / Adding empty rows)
%!test
%! g  = [1; 2];
%! gh = [1; 2];
%! [C, order] = confusionmat (g, gh, "Order", [1; 2; 4]);
%! assert (C, [1 0 0; 0 1 0; 0 0 0]);
%! assert (order, [1; 2; 4]);
## Test 23: All Mismatch
%!test
%! g  = [1; 1; 1];
%! gh = [2; 2; 2];
%! [C, order] = confusionmat (g, gh);
%! assert (C, [0 3; 0 0]);
%! assert (order, [1; 2]);
## Test 24: Single Class Present
%!test
%! g  = [1; 1; 1];
%! gh = [1; 1; 1];
%! [C, order] = confusionmat (g, gh);
%! assert (C, 3);
%! assert (order, 1);
%!error
%! confusionmat ([1; 2], {'A'; 'B'})
%!error
%! confusionmat ('A', [1])
%!error
%! confusionmat ([1; 2; 3], [1; 2])
%!error
%! confusionmat ([1; 2], [1; 2], "Order", {'A'; 'B'})
%!error
%! confusionmat ({'A'}, {'A'}, "Order", [1])
%!error
%! confusionmat (eye(2), eye(2))
