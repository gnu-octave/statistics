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
## @deftypefn {Private Function} [@var{p_net}, @var{p_val}] = exact2xkCT (@var{ct}, @var{weights}, @var{rsstat})
##
## Compute the exact p-value for a 2-by-K contingency table based on the
## network algorithm.
##
## Reference:  Cyrus R. Mehta & Nitin R. Patel (1980) A network algorithm for
## the exact treatment of the 2×k contingency table, Communications in
## Statistics - Simulation and Computation, 9:6, 649-664,
## DOI: 10.1080/03610918008812182
##
## @end deftypefn

function [p_net, p_val] = exact2xkCT (ct, weights, rsstat)

  ## Calculate nodes and arcs
  [nodes, arcs] = build_nodes (ct,weights);

  ## Apply backward induction to nodes
  nodes = backward_induce (nodes,arcs);

  ## Forward scan the network to get p-values
  p_val = forward_scan (nodes, arcs, rsstat);

  ## Calculate p-values
  TP = nodes{4,1};
  p_val = p_val / TP;
  p_net = p_val(2) + min(p_val(1), p_val(3));
endfunction

## Calculate structures describing nodes and arcs
function [nodes, arcs] = build_nodes (ct, weights)
  column = size (ct, 2);      ## number of columns in contigency table
  rowsum = sum (ct, 2);       ## sum of rows
  colsum = sum (ct, 1);       ## sum of columns
  oldnodes = zeros(1,2);      ## nodes added during last pass
  oldlo = 0;                  ## min possible sum so far
  oldhi = 0;                  ## max possible sum so far
  oldnn = 1;                  ## node numbers (row numbers) from last pass
  ctsum = rowsum(1);          ## sum of entries in first row
  nodecount = 1;              ## current node count
  ## Initialize cell structures for nodes and arcs
  nodes = cell(4, column+1);  ## to hold nodes
  nodes{1,1} = zeros (1,2);   ## n-by-2 array, n = # of nodes, row = [j,mj]
  nodes{2,column+1} = 0;      ## n-vector of longest path to end from here
  nodes{3,column+1} = 0;      ## n-vector of shortest path to end from here
  nodes{4,column+1} = 1;      ## n-vector of total probability to end from here
  arcs = cell(3, column);     ## to hold arcs
  ##      row 1:  n-by-2 array, n = # of connections, row = pair connected
  ##      row 2:  n-vector of arc lengths
  ##      row 3:  n-vector of arc probabilities
  for j = 1:column
    ## Find nodes possible at the next step
    nj = colsum(j);
    lo = max (oldlo, ctsum - sum (colsum(j+1:end)));
    hi = min (ctsum, oldhi + nj);
    newnodes = zeros (hi - lo + 1,2);
    newnodes(:,1) = j;
    newnodes(:,2) = (lo:hi)';
    newnn = 1:size (newnodes,1);
    nodecount = nodecount + size (newnodes, 1);
    nodes{1,j+1} = newnodes;
    ## Find arcs possible to the next step
    [a0, a1] = meshgrid (oldnn, newnn);
    a0 = a0(:);
    a1 = a1(:);
    oldsum = oldnodes(a0,2);
    newsum = newnodes(a1,2);
    xj = newsum - oldsum;
    ok = (xj >= 0) & (xj <= nj);
    arcs{1,j} = [a0(ok) a1(ok)];  ## arc connections
    xj = xj(ok);
    arcs{2,j} = weights(j) * xj;
    pj = exp (gammaln (nj + 1) - gammaln (xj + 1) - gammaln (nj - xj + 1));
    arcs{3,j} = pj;               ## arc probabilities
    ## Update data structures
    oldlo = lo;
    oldhi = hi;
    oldnodes = newnodes;
    oldnn = newnn;
  endfor
endfunction

## Calculate backward induction by adding information to NODES array
function nodes = backward_induce (nodes, arcs)
  ## initialize for final node
  column = size (nodes,2) - 1;
  startSP = zeros (1);
  startLP = startSP;
  startTP = ones (1);
  for j = column:-1:1
    ## destination nodes are previous start nodes
    endSP = startSP;
    endLP = startLP;
    endTP = startTP;
    ## get new start nodes and information about them
    a = arcs{1,j};
    startmax = max(a(:,1));
    startSP = zeros(startmax,1);
    startLP = startSP;
    startTP = startSP;
    arclen = arcs{2,j};
    arcprob = arcs{3,j};
    for nodenum = 1:startmax
      % for each start node, compute SP, LP, TP
      k1 = find(a(:,1) == nodenum);
      k2 = a(k1,2);
      startLP(nodenum) = max(arclen(k1) + endLP(k2));
      startSP(nodenum) = min(arclen(k1) + endSP(k2));
      startTP(nodenum) = sum(arcprob(k1) .* endTP(k2));
    endfor
    ## store information about nodes at this level
    nodes{2,j} = startLP;
    nodes{3,j} = startSP;
    nodes{4,j} = startTP;
  endfor
endfunction

## Get p-values by forward scanning the network
function p_val = forward_scan (nodes, arcs, rsstat)
  NROWS = 50;
  p_val = zeros(3,1);    ## [Prob<T, Prob=T, Prob>T]
  stack = zeros(NROWS, 4);
  stack(:,1) = Inf;
  stack(1,1) = 1;        ## level of current node
  stack(1,2) = 1;        ## number at this level of current node
  stack(1,3) = 0;        ## length so far to this node
  stack(1,4) = 1;        ## probability so far of reaching this node
  N = size (stack, 1);
  i1 = 0; i2 = 0; i3 = 0;
  while (1)
    ## Get next lowest level node to process
    minlevel = min(stack((stack(1:N)>0)));
    if (isinf (minlevel))
      break;
    endif
    sp = find (stack(1:N) == minlevel);
    sp = sp(1);
    L = stack(sp,1);
    J = stack(sp,2);
    pastL = stack(sp,3);
    pastP = stack(sp,4);
    stack(sp,1) = Inf;
    ## Get info for arcs at level L and their target nodes
    LP = nodes{2,L+1};
    SP = nodes{3,L+1};
    TP = nodes{4,L+1};
    aj = arcs{1,L};
    arclen = arcs{2,L};
    arcprob = arcs{3,L};
    ## Look only at arcs from node J
    seps = sqrt (eps);
    arows = find (aj(:,1) == J)';
    for k = arows
      tonode = aj(k,2);
      thisL = arclen(k);
      thisP = pastP * arcprob(k);
      len = pastL + thisL;
      ## No paths from node J are significant
      if (len + LP(tonode) < rsstat - seps)
        p_val(1) = p_val(1) + thisP * TP(tonode);
      ## All paths from node J are significant
      elseif (len + SP(tonode) > rsstat + seps)
        p_val(3) = p_val(3) + thisP * TP(tonode);
      ## Single match from node J
      elseif (SP(tonode) == LP(tonode))
        p_val(2) = p_val(2) + thisP * TP(tonode);
      ## Match node J with another already stored node
      else
        ## Find a stored node that matches this one
        r = find(stack(:,1) == L+1);
        if (any (r))
          r = r(stack(r,2) == tonode);
          if (any (r))
             r = r(abs (stack(r,3) - len) < seps);
          endif
        endif
        ## If any one is found, merge node J with it
        if (any (r))
          sp = r(1);
          stack(sp,4) = stack(sp,4) + thisP;
          i1 = i1 + 1;
        ## Otherwise add a new node
        else
          z = find(isinf(stack(:,1)));
          if (isempty (z))
             i2 = i2 +1;
             block = zeros (NROWS, 4);
             block(:,1) = Inf;
             stack = [stack; block];
             sp = N + 1;
             N = N + NROWS;
          else
             i3 = i3 + 1;
             sp = z(1);
          endif
          stack(sp,1) = L + 1;
          stack(sp,2) = tonode;
          stack(sp,3) = len;
          stack(sp,4) = thisP;
        endif
      endif
    endfor
  endwhile
endfunction
