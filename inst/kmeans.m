## Copyright (C) 2011 Soren Hauberg <soren@hauberg.org>
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

function [classes, centers] = kmeans (data, k, varargin)
  ## Input checking
  if (!ismatrix (data) || !isreal (data))
    error ("kmeans: first input argument must be a DxN real data matrix");
  endif
  if (!isscalar (k))
    error ("kmeans: second input argument must be a scalar");
  endif
  
  [N, D] = size (data);
  
  ## (so far) Harcoded options
  maxiter = Inf;
  start = "sample";
  
  ## Find initial clusters
  switch (lower (start))
    case "sample"
      idx = randperm (N) (1:k);
      centers = data (idx, :);
    otherwise
      error ("kmeans: unsupported initial clustering parameter");
  endswitch
  
  ## Run the algorithm
  D = zeros (N, k);
  iterations = 0;
  prevcenters = centers;
  while (true)
    ## Compute distances
    for i = 1:k
      D (:, i) = sum (( data - repmat (centers (i, :), N, 1)).^2, 2);
    endfor
    
    ## Classify
    [tmp, classes] = min (D, [], 2);
    
    ## Recompute centers
    for i = 1:k
      centers (i, :) = mean (data (classes == i, :));
    endfor
    
    ## Check for convergence
    iterations++;
    if (all (centers (:) == prevcenters (:)) || iterations >= maxiter)
      break;
    endif
    prevcenters = centers;
  endwhile
endfunction

%!demo
%! ## Generate a two-cluster problem
%! C1 = randn (100, 2) + 1;
%! C2 = randn (100, 2) - 1;
%! data = [C1; C2];
%!
%! ## Perform clustering
%! [idx, centers] = kmeans (data, 2);
%!
%! ## Plot the result
%! figure
%! plot (data (idx==1, 1), data (idx==1, 2), 'ro');
%! hold on
%! plot (data (idx==2, 1), data (idx==2, 2), 'bs');
%! plot (centers (:, 1), centers (:, 2), 'kv', 'markersize', 10);
%! hold off

