## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the Free
## Software Foundation; either version 3 of the License, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {statistics} {@var{node} =} __build_kdtree__ (@var{indices}, @var{depth}, @var{X}, @var{bucket_size})
##
## Build a Kd-tree over the rows of @var{X} named by @var{indices}.  Shared by
## @code{knnsearch}, @code{rangesearch} and @code{KDTreeSearcher}.
##
## A node is a leaf, carrying an @qcode{indices} field, or a split, carrying
## @qcode{axis}, @qcode{split_value}, @qcode{left} and @qcode{right}.  The
## split axis cycles with the depth.
##
## @end deftypefn

function node = __build_kdtree__ (indices, depth, X, bucket_size)

  if (length (indices) <= bucket_size)
    node = struct ('indices', indices);
  else
    k = size (X, 2);
    axis = mod (depth, k) + 1;
    values = X(indices, axis);
    sorted_values = sort (values);
    median_idx = floor ((length (indices) + 1) / 2);
    split_value = sorted_values(median_idx);
    left_indices = indices(values <= split_value);
    right_indices = indices(values > split_value);

    ## Points sharing the split value all go left, so a set that is constant
    ## on this axis cannot be divided here.  Without this the recursion hands
    ## itself the same set forever, which more than BUCKET_SIZE identical rows
    ## used to do.
    if (isempty (left_indices) || isempty (right_indices))
      node = struct ('indices', indices);
      return;
    endif

    left_node = __build_kdtree__ (left_indices, depth + 1, X, bucket_size);
    right_node = __build_kdtree__ (right_indices, depth + 1, X, bucket_size);
    node = struct ('axis', axis, 'split_value', split_value, ...
                   'left', left_node, 'right', right_node);
  endif

endfunction
