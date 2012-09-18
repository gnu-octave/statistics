%% Copyright (c) 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
%%
%%    This program is free software: you can redistribute it and/or modify
%%    it under the terms of the GNU General Public License as published by
%%    the Free Software Foundation, either version 3 of the License, or
%%    any later version.
%%
%%    This program is distributed in the hope that it will be useful,
%%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%    GNU General Public License for more details.
%%
%%    You should have received a copy of the GNU General Public License
%%    along with this program. If not, see <http://www.gnu.org/licenses/>.

%% -*- texinfo -*-
%% @deftypefn {Function File} {@var{h} = } dendogram (@var{tree})
%% Plots a dendogram using the output of function @command{linkage}.
%%
%% TODO: This function needs your help to be finished! Leafs must be ordered to
%%       prevent corssing of edges.
%% @seealso{linkage}
%% @end deftypefn

## TODO
# Order leafs nodes to avoid edge crossing.

function h = dendogram (tree)

  [m d] = size (tree);
  if d != 3
    error ("Input data must be a tree as returned by function linkage.")
  end
  n = m + 1;

  nc = max(tree(:,1:2)(:));

  % Vector with the horizontal and vertical position of each cluster
  p = zeros (nc,2);

##### This is an ugly hack and is a partial solution. If you know how to oder
##### the leafs do it. Please do not copy from dendogram.m of Mathworks!

  %% Order the leafs
  % find which clusters contain th eleafs
  tf = ismember (tree(:,1:2),1:n);
  [idx_i,idx_j] = find (tf);
  ind0 = find(tf(:));

  % Assign a position between 1:n to the leafs according to clustering
  % I am trying to avoid crossings. Better way? sure!.
  i = 1;
  j = 0;
  while j < n

    % If the current leaf hasn't a position assigned, assign one.
    if p(tree(idx_i(i),idx_j(i)),1) == 0
     j+=1;
     p(tree(idx_i(i),idx_j(i)),1) = j;
    end

    % Check if the current cluster has another leaf. If it has and the leaf
    % doesn't have a position, assign one.
    next = tree(idx_i(i),setdiff(1:2,idx_j(i)));
    if next <= n && p(next,1) == 0
      j+=1;
      p(next,1) = j;
    end

    % Find the id of the current cluster and check if any other cluster
    % that will be merged with this one has a leaf.
    c_id = idx_i(i) + n;
    [idx,~] = find (tree(:,1:2)==c_id);
    idx2 = find (tree(idx,1:2)<=n);
    if ~isempty (idx2) && p(tree(idx,idx2),1) == 0
      j+=1;
      p(tree(idx,idx2),1) = j;
    end
    i+=1;
  end

##### End of the hack

  % Compute the horizontal position, begin-end
  % and vertical position of all clusters.
  for i = 1:m
    p(n+i,1)   = mean (p(tree(i,1:2),1));
    p(n+i,2)   = tree(i,3);
    x(i,1:2) = p(tree(i,1:2),1);
  end



  clf
  % plot horizontal lines
  tmp = line (x', tree(:,[3 3])');

  % plot vertical lines
  for i=1:nc
    [ind,~] = find (tree(:,1:2)==i);
    tmp = line (p([i; i],1),[p(i,2); tree(ind,3)]);
  end

  xticks = 1:n;
  xl_txt = arrayfun (@num2str, tree(:,1:2)(ind0),"uniformoutput",false);
  set (gca,"xticklabel",xl_txt,"xtick",xticks);
  axis ([0.5 n+0.5 0 max(tree(:,3))+0.1*min(tree(:,3))]);

endfunction

%!demo
%! y      = [4 5; 2 6; 3 7; 8 9; 1 10];
%! y(:,3) = 1:5;
%! dendogram(y);
