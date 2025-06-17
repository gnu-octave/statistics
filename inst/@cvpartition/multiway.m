## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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

function [groupindex, partition, groupsizes] = multiway (numbers, num_parts, varargin)
  if nargin < 2
    num_parts = 2;
  endif

  method = "completeKK";
  if nargin > 2
    if mod (length (varargin), 2) ~= 0
      error (strcat ("Invalid number of arguments. Additional arguments", ...
                     " must be name-value pairs."));
    endif
    for i = 1:2:length (varargin)
      if strcmpi (varargin{i}, "method")
        method = varargin{i+1};
      else
        error ("Unknown parameter: %s", varargin{i});
      endif
    endfor
  endif

  switch lower(method)
    case "completekk"
      [groupindex, partition, groupsizes] = complete_karmarkar_karp (numbers, ...
                                                                     num_parts);
    otherwise
      error ("multiway: unsupported method '%s'", method);
  endswitch
endfunction

function [groupindex, partition, groupsizes] = complete_karmarkar_karp (numbers, num_parts)
  if isempty (numbers)
    partition = cell (1, num_parts);
    for i = 1:num_parts
      partition{i} = [];
    endfor
    groupindex = zeros (size (numbers));
    groupsizes = zeros (1, num_parts);
    return;
  endif

  initial_nodes = cell (1, numel (numbers));
  for i = 1:numel (numbers)
    part = cell (1, num_parts);
    part_sums = zeros (1, num_parts);
    part{end} = numbers(i);
    part_sums(end) = numbers(i);
    initial_nodes{i} = struct ( ...
      'parts', {part}, ...
      'sums', part_sums, ...
      'diff', numbers(i), ...
      'id', i ...
    );
  endfor

  processing_stack = {initial_nodes};
  best_diff = Inf;
  best_node = [];
  id_counter = numel (numbers);

  while ~isempty (processing_stack)
    current_nodes = processing_stack{end};
    processing_stack(end) = [];

    min_possible_diff = calculate_lower_bound (current_nodes, num_parts);
    if min_possible_diff >= best_diff
      continue;
    endif

    if numel (current_nodes) == 1
      current_diff = current_nodes{1}.diff;
      if current_diff < best_diff
        best_diff = current_diff;
        best_node = current_nodes{1};
        if current_diff == 0
          break;
        endif
      endif
      continue;
    endif

    [~, order] = sort (cellfun (@(x) x.diff, current_nodes));
    current_nodes = current_nodes(order);
    node1 = current_nodes{1};
    node2 = current_nodes{2};
    remaining_nodes = current_nodes(3:end);

    new_nodes = combine_partitions (node1, node2, id_counter);
    id_counter = id_counter + numel (new_nodes);

    candidate_sets = cell (1, numel (new_nodes));
    for i = 1:numel (new_nodes)
      candidate_sets{i} = [remaining_nodes, new_nodes{i}];
    endfor

    candidate_bounds = cellfun (@(x) calculate_lower_bound (x, num_parts), ...
                                candidate_sets);
    [~, order] = sort (candidate_bounds, 'ascend');
    processing_stack = [processing_stack, candidate_sets(order)];
  endwhile

  partition = best_node.parts;
  groupsizes = best_node.sums;

  if nargout >= 1
    idx_cell = convert_to_indices (partition, numbers);
    groupindex = zeros (size (numbers));
    for j = 1:num_parts
      groupindex(idx_cell{j}) = j;
    endfor
  endif
endfunction

function indices = convert_to_indices (parts, original_numbers)
  indices = cell (size (parts));
  numbers_copy = original_numbers(:)';

  for i = 1:numel (parts)
    idxs = [];
    for val = parts{i}
      pos = find (numbers_copy == val, 1);
      if isempty (pos)
        error ("Value %d not found during index conversion", val);
      endif
      idxs = [idxs, pos];
      numbers_copy(pos) = NaN;
    endfor
    indices{i} = sort (idxs);
  endfor
endfunction

function new_nodes = combine_partitions (node1, node2, start_id)
  parts1 = node1.parts;
  parts2 = node2.parts;
  num_parts = numel (parts1);

  all_perm = perms (1:num_parts);
  num_perm = size (all_perm, 1);
  new_nodes = cell (1, num_perm);

  for idx = 1:num_perm
    p_row = all_perm(idx, :);
    combined_parts = cell (1, num_parts);
    combined_sums = zeros (1, num_parts);

    for i = 1:num_parts
      combined_parts{i} = [parts1{p_row(i)}, parts2{i}];
      combined_sums(i) = sum (combined_parts{i});
    endfor

    part_diff = max (combined_sums) - min (combined_sums);

    new_node = struct ( ...
      'parts', {combined_parts}, ...
      'sums', combined_sums, ...
      'diff', part_diff, ...
      'id', start_id + idx - 1 ...
    );

    new_nodes{idx} = new_node;
  endfor
endfunction

function bound = calculate_lower_bound (nodes, num_parts)
  all_sums = cellfun (@(x) x.sums, nodes, "UniformOutput", false);
  all_sums = [all_sums{:}];
  total = sum (all_sums);
  max_val = max (all_sums);
  bound = max_val - (total - max_val) / (num_parts - 1);
endfunction

%!test
tic;
numbers = [4, 5, 6, 7, 8];
num_parts = 2;
[groupindex, partition, groupsizes] = multiway (numbers, num_parts, "method", "completeKK");
assert (isequal (sort (cellfun (@sum, partition)), sort ([15, 15])))
toc

%! test
tic;
numbers = [1, 2, 3, 4, 5, 6];
num_parts = 3;
[groupindex, partition, groupsizes] = multiway (numbers, num_parts, "method", "completeKK");
assert (isequal (sort (cellfun (@sum, partition)), sort ([7, 7, 7])))
toc

%! test
tic;
numbers = [24, 21, 18, 17, 12, 11, 8, 2];
num_parts = 3;
[groupindex, partition, groupsizes] = multiway (numbers, num_parts, "method", "completeKK");
assert (isequal (sort (cellfun (@sum, partition)), sort ([38, 38, 37])))
toc

%! test
tic;
numbers = [127, 125, 122, 105, 87, 75, 68, 64, 30, 22];
num_parts = 4;
[groupindex, partition, groupsizes] = multiway (numbers, num_parts, "method", "completeKK");
assert (isequal (sort (cellfun (@sum, partition)), sort ([211, 202, 209, 203])))
toc