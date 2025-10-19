## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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
## @deftypefn  {statistics} {@var{groupindex} =} multiway (@var{numbers}, @var{num_parts})
## @deftypefnx {statistics} {@var{groupindex} =} multiway (@var{numbers}, @var{num_parts}, @var{method})
## @deftypefnx {statistics} {[@var{groupindex}, @var{partition}] =} multiway (@dots{})
## @deftypefnx {statistics} {[@var{groupindex}, @var{partition}, @var{groupsizes}] =} multiway (@dots{})
##
## Solve the multiway number partitioning problem.
##
## @code{@var{groupindex} = multiway (@var{numbers}, @var{num_parts})} splits
## a set of numbers in @var{numbers} into a number of subsets specified in
## @var{num_parts} such that the sums of the subsets are nearly as equal as
## possible and returns a vector of group indices in @var{groupindex} with each
## index corresponding to the set of numbers provided as input.
##
## @itemize
## @item @var{numbers} is a vector of positive real numbers to be partitioned.
## @item @var{num_parts} is a positive integer scalar specifying the number of
## partitions (subsets) to split the numbers into.
## @end itemize
##
## @code{@var{groupindex} = multiway (@var{numbers}, @var{num_parts},
## @var{method})} also specifies the algorithm used for partitioning the set of
## numbers.  By default, @code{multiway} uses the complete Karmarkar-Karp
## algorithm, when the set of numbers contains up to 10 elements and the
## requested number of subsets does not exceed 5, otherwise it defaults to the
## greedy algorithm, which is optimized for speed, but may not return the
## optimal partitioning.  The following methods are supported:
##
## @itemize
## @item @qcode{'greedy'} (Greedy algorithm)
## @item @qcode{'completeKK'} (Complete Karmarkar-Karp algorithm)
## @end itemize
##
## The @code{multiway} function may return up to three output arguments
## described below:
##
## @itemize
## @item
## @var{groupindex}: A vector of the same length as @var{numbers} containing
## the group index (from 1 to @var{num_parts}) for each number.
## @item
## @var{partition}: A cell array of length @var{num_parts} with each cell
## containing the numbers assigned to that partition.
## @item
## @var{groupsizes}: A vector of the sums of the numbers in each partition.
## @end itemize
##
## Example:
## @example
## @group
## numbers = [4, 5, 6, 7, 8];
## num_parts = 2;
## [groupindex, partition, groupsizes] = multiway (numbers, num_parts);
## @end group
## @end example
##
## @seealso{cvpartition}
## @end deftypefn

function [gindex, partition, gsize] = multiway (numbers, num_parts, method)

  if (nargin < 2)
    error ("multiway: too few input arguments.");
  elseif (nargin == 2)
    if (numel (numbers) <= 10 && num_parts <= 5)
      method = 'completeKK';
    else
      method = 'greedy';
    endif
  elseif (! (ischar (method) && isvector (method)))
    error ("multiway: METHOD value must be a character vector.");
  endif

  ## Validate numbers vector
  if (! isvector (numbers) || isempty (numbers))
    error ("multiway: NUMBERS must be a non-empty vector.");
  endif
  if (! isnumeric (numbers) || any (isnan (numbers)))
    error ("multiway: NUMBERS must be numeric and cannot contain NaN values.");
  endif
  if (any (numbers < 0))
    error ("multiway: NUMBERS must be non-negative.");
  endif

  ## Validate number of partitions
  if (! isscalar (num_parts))
    error ("multiway: NUM_PARTS must be a scalar.");
  endif
  if (! isnumeric (num_parts) || ! isreal (num_parts))
    error ("multiway: NUM_PARTS must be a real numeric value.");
  endif
  if (num_parts < 1 || fix (num_parts) != num_parts)
    error ("multiway: NUM_PARTS must be a positive integer.");
  endif
  if (num_parts > numel (numbers))
    error (strcat ("multiway: NUM_PARTS cannot be greater than", ...
                   " number of elements in NUMBERS."));
  endif

  ## Select method
  switch (lower (method))
    case 'completekk'
      [gindex, partition, gsize] = completeKK (numbers, num_parts);
    case 'greedy'
      [gindex, partition, gsize] = greedy (numbers, num_parts);
    otherwise
      error ("multiway: unsupported method '%s'.", method);
  endswitch
endfunction

function [groupindex, partition, groupsizes] = greedy (numbers, num_parts)
  [sorted_numbers, sorted_indices] = sort (numbers, 'descend');
  n = numel (sorted_numbers);

  partition = cell (1, num_parts);
  sums = zeros (1, num_parts);

  group_assignment = zeros (1, n);
  for i = 1:n
    [min_sum, min_idx] = min (sums);
    partition{min_idx}(end+1) = sorted_numbers(i);
    sums(min_idx) = min_sum + sorted_numbers(i);
    group_assignment(i) = min_idx;
  endfor

  groupindex = zeros (size (numbers));
  groupindex(sorted_indices) = group_assignment;

  groupsizes = sums;
  if (iscolumn (numbers))
    groupsizes = groupsizes';
  endif
endfunction

function [groupindex, partition, groupsizes] = completeKK (numbers, num_parts)
  [gidx_g, part_g, gsize_g] = greedy (numbers, num_parts);
  best_diff = max (gsize_g) - min (gsize_g);
  best_partition = part_g;
  groupindex = gidx_g;

  if (best_diff == 0)
    partition = best_partition;
    groupsizes = gsize_g;
    if (iscolumn (numbers))
      groupsizes = groupsizes';
    endif
    return;
  endif

  sorted_numbers = numbers(:).';
  [sorted_numbers, sorted_indices] = sort (sorted_numbers, 'descend');
  n = numel (sorted_numbers);
  total = sum (sorted_numbers);
  average = total / num_parts;

  function rec (k, current_sums, current_assign)
    if (k > n)
      this_diff = max (current_sums) - min (current_sums);
      if (this_diff < best_diff)
        best_diff = this_diff;
        for jj = 1:num_parts
          best_partition{jj} = sorted_numbers(current_assign == jj);
        endfor
      endif
      return;
    endif

    rem_max = sorted_numbers(k);
    this_min = min (current_sums);
    lb_max = max (max (current_sums), this_min + rem_max);
    lb_diff = max (0, lb_max - average);
    if (lb_diff >= best_diff)
      return;
    endif

    [~, perm] = sort (current_sums);

    for jj = 1:num_parts
      j = perm (jj);
      new_sums = current_sums;
      new_sums (j) += sorted_numbers (k);

      new_max = max (new_sums);
      if (k < n)
        next_rem_max = sorted_numbers (k + 1);
        sub_min = min (new_sums);
        sub_lb_max = max (new_max, sub_min + next_rem_max);
      else
        sub_lb_max = new_max;
      endif
      sub_lb_diff = max (0, sub_lb_max - average);
      if (sub_lb_diff >= best_diff)
        continue;
      endif

      new_assign = current_assign;
      new_assign (k) = j;
      rec (k + 1, new_sums, new_assign);
    endfor
  endfunction

  rec (1, zeros (1, num_parts), zeros (1, n));

  partition = best_partition;
  groupsizes = cellfun (@sum, partition);

  idx_cell = convert_to_indices (partition, numbers);
  groupindex = zeros (size (numbers));
  for j = 1:num_parts
    groupindex(idx_cell{j}) = j;
  endfor

  if (iscolumn (numbers))
    groupsizes = groupsizes';
  endif
endfunction

function indices = convert_to_indices (parts, original_numbers)
  indices = cell (size (parts));
  numbers_copy = original_numbers(:)';

  for i = 1:numel (parts)
    idxs = [];
    for val = parts{i}
      pos = find (numbers_copy == val, 1);
      if (isempty (pos))
        error ("Value %g not found during index conversion", val);
      endif
      idxs = [idxs, pos];
      numbers_copy(pos) = NaN;
    endfor
    indices{i} = idxs;
  endfor
endfunction

## Test completeKK method
%!test
%! numbers = [4, 5, 6, 7, 8];
%! num_parts = 2;
%! [groupindex, partition, groupsizes] = multiway (numbers, num_parts, "completeKK");
%! assert (sort (cellfun (@sum, partition)), sort ([15, 15]));
%!test
%! numbers = [1, 2, 3, 4, 5, 6];
%! num_parts = 3;
%! [groupindex, partition, groupsizes] = multiway (numbers, num_parts, "completeKK");
%! assert (sort (cellfun (@sum, partition)), sort ([7, 7, 7]));
%!test
%! numbers = [24, 21, 18, 17, 12, 11, 8, 2];
%! num_parts = 3;
%! [groupindex, partition, groupsizes] = multiway (numbers, num_parts, "completeKK");
%! assert (sort (cellfun (@sum, partition)), sort ([38, 38, 37]));
%!test
%! numbers = [10, 10, 10];
%! num_parts = 3;
%! [~, partition] = multiway (numbers, num_parts, "completeKK");
%! assert (sort (cellfun (@sum, partition)), [10, 10, 10]);
%!test
%! numbers = 1:10;
%! num_parts = 2;
%! [~, partition] = multiway (numbers, num_parts, "completeKK");
%! assert (sort (cellfun (@sum, partition)), [27, 28]);

## Test greedy method
%!test
%! numbers = [4, 5, 6, 7, 8];
%! num_parts = 2;
%! [groupindex, partition, groupsizes] = multiway (numbers, num_parts, "greedy");
%! assert (sort (cellfun (@sum, partition)), sort ([13, 17]));
%!test
%! numbers = [1, 2, 3, 4, 5, 6];
%! num_parts = 3;
%! [groupindex, partition, groupsizes] = multiway (numbers, num_parts, "greedy");
%! assert (sort (cellfun (@sum, partition)), sort ([7, 7, 7]));
%!test
%! numbers = [10, 7, 5, 5, 6, 4, 10, 11, 12, 9, 10, 4, 3, 4, 5];
%! num_parts = 4;
%! [groupindex, partition, groupsizes] = multiway (numbers, num_parts, "greedy");
%! assert (sort (cellfun (@sum, partition)), sort ([27, 27, 27, 24]));
%!test
%! numbers = [24, 21, 18, 17, 12, 11, 8, 2];
%! num_parts = 3;
%! [groupindex, partition, groupsizes] = multiway (numbers, num_parts, "greedy");
%! assert (sort (cellfun (@sum, partition)), sort ([35, 37, 41]));
%!test
%! numbers = [10, 10, 10];
%! num_parts = 3;
%! [~, partition] = multiway (numbers, num_parts, "greedy");
%! assert (sort (cellfun (@sum, partition)), [10, 10, 10]);
%!test
%! numbers = 1:10;
%! num_parts = 2;
%! [~, partition] = multiway (numbers, num_parts, "greedy");
%! assert (sort (cellfun (@sum, partition)), [27, 28]);

## Test algorithm switch
%!test
%! grpidx_ckk = multiway ([3 2 4 3 9 3 64], 3);
%! grpidx_greedy = multiway ([3 2 4 3 9 3 64], 3, 'greedy');
%! assert (isequal (grpidx_ckk, grpidx_greedy), false);

## Test column vector input
%!test
%! numbers = [4; 5; 6; 7; 8];
%! num_parts = 2;
%! [groupindex, partition, groupsizes] = multiway (numbers, num_parts, "completeKK");
%! assert (iscolumn (groupindex), true)
%! assert (iscolumn (groupsizes), true);
%! assert (sort (cellfun (@sum, partition)), sort ([15, 15]));
%! numbers = [4; 5; 6; 7; 8];
%! num_parts = 2;
%! [groupindex, partition, groupsizes] = multiway (numbers, num_parts, "greedy");
%! assert (iscolumn (groupindex), true)
%! assert (iscolumn (groupsizes), true);
%! assert (sort (cellfun (@sum, partition)), sort ([13, 17]));

## Test row vector input
%!test
%! numbers = [4, 5, 6, 7, 8];
%! num_parts = 2;
%! [groupindex, partition, groupsizes] = multiway (numbers, num_parts, "completeKK");
%! assert (isrow (groupindex), true)
%! assert (isrow (groupsizes), true);
%! assert (sort (cellfun (@sum, partition)), sort ([15, 15]));
%!test
%! numbers = [4, 5, 6, 7, 8];
%! num_parts = 2;
%! [groupindex, partition, groupsizes] = multiway (numbers, num_parts, "greedy");
%! assert (isrow (groupindex), true)
%! assert (isrow (groupsizes), true);
%! assert (sort (cellfun (@sum, partition)), sort ([13, 17]));

## Test input validation
%!error<multiway: too few input arguments.> multiway ()
%!error<multiway: too few input arguments.> multiway ([1, 2])
%!error<multiway: METHOD value must be a character vector.> ...
%! multiway ([1, 2, 3], 2, 1)
%!error<multiway: NUMBERS must be a non-empty vector.> multiway ([], 2)
%!error<multiway: NUMBERS must be a non-empty vector.> multiway (ones (2, 2), 2)
%!error<multiway: NUMBERS must be numeric and cannot contain NaN values.> ...
%! multiway ({1, 2, 3}, 2)
%!error<multiway: NUMBERS must be non-negative.> multiway ([1, -2, 3], 2)
%!error<multiway: NUMBERS must be numeric and cannot contain NaN values.> ...
%! multiway ([1, 2, NaN], 2)
%!error<multiway: NUM_PARTS must be a scalar.> multiway ([1,2,3], [1,2])
%!error<multiway: NUM_PARTS must be a real numeric value.> ...
%! multiway ([1, 2, 3], "2")
%!error<multiway: NUM_PARTS must be a positive integer.> multiway ([1, 2, 3], 0)
%!error<multiway: NUM_PARTS must be a positive integer.> multiway ([1, 2, 3], 1.5)
%!error<multiway: NUM_PARTS must be a positive integer.> multiway ([1, 2, 3], -1)
%!error<multiway: NUM_PARTS cannot be greater than number of elements in NUMBERS.> ...
%! multiway ([1, 2], 3)
%!error <multiway: unsupported method 'greedyalgo'.> ...
%! multiway ([1,2,3], 2, "greedyalgo")
