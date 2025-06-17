function [groupindex, partition, groupsizes] = complete_karmarkar_karp (numbers, num_parts, varargin)
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

  if ~strcmpi (method, "completeKK")
    error ("Only completeKK method is currently supported");
  endif

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

    candidate_bounds = cellfun (@(x) calculate_lower_bound (x, num_parts), candidate_sets);
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

function part_struct = format_partition (node, original_numbers, return_indices)
  parts = node.parts;
  sums = node.sums;

  if return_indices
    parts = convert_to_indices (parts, original_numbers);
  endif

  part_struct = struct ( ...
    'partition', {parts}, ...
    'sums', sums, ...
    'difference', max (sums) - min (sums) ...
  );
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
  new_nodes = {};

  all_perm = perms (1:num_parts);
  unique_perm = unique (all_perm, 'rows');

  for p = unique_perm'
    combined_parts = cell (1, num_parts);
    combined_sums = zeros (1, num_parts);

    for i = 1:num_parts
      combined_parts{i} = sort ([parts1{p(i)}, parts2{i}]);
      combined_sums(i) = sum (combined_parts{i});
    endfor

    part_diff = max (combined_sums) - min (combined_sums);

    new_node = struct ( ...
      'parts', {combined_parts}, ...
      'sums', combined_sums, ...
      'diff', part_diff, ...
      'id', start_id ...
    );

    new_nodes{end + 1} = new_node;
    start_id = start_id + 1;
  endfor
endfunction

function bound = calculate_lower_bound (nodes, num_parts)
  all_sums = cellfun (@(x) x.sums, nodes, 'UniformOutput', false);
  all_sums = [all_sums{:}];
  total = sum (all_sums);
  max_val = max (all_sums);
  bound = max_val - (total - max_val) / (num_parts - 1);
endfunction