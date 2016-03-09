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
%% @deftypefn {Function File} {@var{p} = } dendrogram (@var{tree})
%% @deftypefnx {Function File} {@var{p, t} = } dendrogram (@var{tree})
%% @deftypefnx {Function File} {@var{p, t, perm} = } dendrogram (@var{tree})
%% Plots a dendrogram using the output of function @command{linkage}.
%%
%% t is a vector containing the leaf node number for each object in the
%% original dataset.  For now, all objects are leaf nodes.
%%
%% perm is the permutation of the input objects used to display the
%% dendrogram, in left-to-right order.
%%
%% TODO: Return handle to lines to set properties
%% TODO: Rescale the plot automatically based on data.
%%
%% @seealso{linkage}
%% @end deftypefn

function [p, t, perm] = dendrogram (tree)

  [m d] = size (tree);
  if d != 3
    error ("Input data must be a tree as returned by function linkage.")
  end
  n = m + 1;

  % t is the leaf node number for all objects in the original dataset.
  % TODO: Add support for collapsing the tree.
  % For now, we always display all objects, so this is the identity map.
  t = (1:m)';

  nc = max(tree(:,1:2)(:));

  % Vector with the horizontal and vertical position of each cluster
  p = zeros (nc,2);

  perm = zeros (n,1);

  %% Ordering by depth-first search
  nodecount = 0;
  nodes_to_visit = nc+1;
  while !isempty(nodes_to_visit)
    currentnode = nodes_to_visit(1);
    nodes_to_visit(1) = [];
    if currentnode > n
      node = currentnode - n;
      nodes_to_visit = [tree(node,[2 1]) nodes_to_visit];
    end

    if currentnode <= n && p(currentnode,1) == 0
      nodecount +=1;
      p(currentnode,1) = nodecount;
      perm(nodecount) = currentnode;
    end

  end

  % Compute the horizontal position, begin-end
  % and vertical position of all clusters.
  for i = 1:m
    p(n+i,1)   = mean (p(tree(i,1:2),1));
    p(n+i,2)   = tree(i,3);
    x(i,1:2) = p(tree(i,1:2),1);
  end

  figure(gcf)
  % plot horizontal lines
  tmp = line (x', tree(:,[3 3])');

  % plot vertical lines
  [~,tf]  = ismember (1:nc, tree(:,1:2));
  [ind,~] = ind2sub (size (tree(:,1:2)), tf);
  y       = [p(1:nc,2) tree(ind,3)];
  tmp     = line ([p(1:nc,1) p(1:nc,1)]',y');

  xticks = 1:n;
  xl_txt = arrayfun (@num2str, perm,"uniformoutput",false);
  set (gca,"xticklabel",xl_txt,"xtick",xticks);
  axis ([0.5 n+0.5 0 max(tree(:,3))+0.1*min(tree(:,3))]);

endfunction

%!demo
%! y      = [4 5; 2 6; 3 7; 8 9; 1 10];
%! y(:,3) = 1:5;
%! figure(gcf); clf;
%! dendrogram(y);

%!demo
%! v = 2*rand(30,1)-1;
%! d = abs(bsxfun(@minus, v(:,1), v(:,1)'));
%! y = linkage (squareform(d,"tovector"));
%! figure(gcf); clf;
%! dendrogram(y);
